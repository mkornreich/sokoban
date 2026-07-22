// Move-OPTIMAL Sokoban solver (Sokoban Variant: keeper must also end on a goal).
//
// Techniques:
// * Compact encoding: boxes = a bitmask over the "live" cells, keeper = a cell
// index (uint8).  Boards with <=32 live cells use a uint32 mask and a packed
// 64-bit lock-free hash ((mask<<8|keeper)<<24 | value). The fast path.  Larger
// boards (up to 128 live cells) switch to a 128-bit mask with a wider lock-free
// hash. The whole search is templated on the mask type so the small case pays
// nothing for the large one.
// * Open-addressed flat hash (linear probing) for best-g / closed, backed by
// 2MB-aligned memory with MADV_HUGEPAGE (the table is hit at random, so TLB
// misses would otherwise dominate) and software prefetch of upcoming slots.
// * Dial bucket queue for the A* open list (integer f).  Fast-path entries pack
// (mask,g,keeper) into one uint64.  Successor h is computed incrementally:
// at expansion h(mask) = f - g exactly, so h(nmask) = f - g - DMIN[from] +
// DMIN[to]. No per-successor loop over boxes.
// * Macro-push moves with move-cost edges: each edge is (shortest keeper walk
// to the push spot) + 1, so g counts real keeper MOVES and the answer is the
// true minimum-move solution, while the state space stays (boxes, keeper).
// * Dead-square table (reverse-reachability from goals) + freeze-deadlock
// pruning. Successors with f >= best-known solution are never queued, and a
// level stops as soon as best <= f (all cheaper levels are already exhausted,
// so best is provably optimal).
// * Workers keep per-thread counters and per-(thread,f) successor chunks that
// are bulk-merged after each wave. No shared-cacheline atomics in the loop.
// * Unsigned throughout: cells, distances, and live-bit indices are uint32 with
// NONE = UINT32_MAX as the single "missing" sentinel (it also compares larger
// than every real distance, so min-updates need no separate presence check).
//
// Reports the optimal move count.  With --solution (-s) it also reconstructs and
// prints an optimal move path (lowercase u/d/l/r = keeper walk, uppercase = push),
// recovered from the completed g-map after the search.
// Usage: ./solver <boardId> [timeoutSec] [--solution|-s] [dbg]
// Build: g++ -O3 -march=native -std=c++17 -pthread solver.cpp -o solver

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <chrono>
#include <thread>
#include <atomic>
#include <algorithm>
#include <mutex>
#include <condition_variable>
#include <functional>
#if defined(__linux__)
#include <sys/mman.h>
#endif
using namespace std;

// ---- boards, rows as digit strings (0 floor,1 wall,2 box,3 keeper,4 goal,5 box+goal,6 keeper+goal)
static unordered_map<string, vector<string>> BOARDS = {
  {"p1", {"111111","103001","102001","110111","100001","104041","111111"}},
  {"p2", {"1111111","1000001","1000001","1002141","1340101","1111111"}},
  {"p3", {"111111111","100010001","100020341","100010001","104010001","111111111"}},
  {"p4", {"1111111","0000014","0000000","0011100","0010000","0201040","0301000"}},
  {"p5", {"111111","110011","100001","142241","100041","113111","111111"}},
  {"p6", {"11111111","10000041","14002231","10010041","11111111"}},
  {"p7", {"1111111111","0011114003","0000010000","0000010010","0010010010","0210000010","0010000014"}},
  {"p8", {"111111","140041","102201","120101","134041","111111"}},
  {"p9", {"111111111","111001111","100000201","101001201","104441301","111111111"}},
  {"p10", {"1111100","1400110","1320011","1102001","0110201","0011001","0001141","0000141","0000141","0000111"}},
  {"p11", {"0010000","0214040","0204000","3211140","0014000"}},
  {"p12", {"11111000","10041000","12101111","14000001","10050501","10501011","11103010","00111110"}},
  {"p13", {"1111111111","1300100441","1020200441","1022211441","1000011441","1111110000"}},
  {"p14", {"0000100000001000","0000100000001000","1111100000011111","0000010000100000","0000001001000000","0000000030000000","0000001001000000","0000010000100000","1111100000011111","0000100000001000","0000100000001000","0000100000401000","0000102000401000","0000102000401000"}},
  {"p15", {"0011111110","1110011110","1002000110","1320200010","1102020010","0110202010","0011024010","0001111010","0000141001","0000144401","0000101401","0000144401","0000111111"}},
};

constexpr uint32_t NONE = UINT32_MAX;   // missing cell / distance / live index
constexpr size_t NPOS = SIZE_MAX;       // hash slot not found

unsigned R, C, N;
vector<uint8_t> WALL, ISGOAL;
vector<uint32_t> NB;            // N*4 neighbor cell or NONE
vector<uint32_t> DMIN;          // min push-distance to any goal, NONE = dead
vector<uint32_t> liveIdx;       // cell -> bit index (0..L-1) or NONE
vector<uint32_t> liveCell;      // bit index -> cell
vector<uint32_t> goalCells;
unsigned L;                     // number of live cells
vector<uint32_t> boxCells0; uint32_t keeper0;   // start boxes / keeper (mask built per width)
const int DR[4] = {-1, 1, 0, 0}, DC[4] = {0, 0, -1, 1};   // signed by nature: direction deltas
const unsigned OPP[4] = {1, 0, 3, 2};

// value layout shared by both hashes: bit 23 = closed, bits 0..22 = best g.
constexpr uint32_t VALMASK = (1u << 24) - 1;
constexpr uint32_t CLOSEDBIT = 1u << 23;
constexpr uint32_t GMASK = (1u << 23) - 1;

inline uint32_t idx(unsigned r, unsigned c) { return r * C + c; }
inline bool inb(int r, int c) { return r >= 0 && r < (int)R && c >= 0 && c < (int)C; }

// 2MB-aligned allocation marked for transparent huge pages: the hash tables are
// gigabytes hit at random, so 4K pages would thrash the TLB.
static void* bigalloc(size_t bytes) {
  const size_t ALIGN = 2ull << 20;
  bytes = (bytes + ALIGN - 1) & ~(ALIGN - 1);
  void* p = aligned_alloc(ALIGN, bytes);
  if (!p) { fprintf(stderr, "out of memory (%zu bytes)\n", bytes); exit(2); }
#if defined(__linux__)
  madvise(p, bytes, MADV_HUGEPAGE);
#endif
  return p;
}

// ---- mask width helpers (overloaded so templated code picks the right builtin) ----
static inline unsigned MPOP(uint32_t m) { return (unsigned)__builtin_popcount(m); }
static inline unsigned MPOP(unsigned long long m) { return (unsigned)__builtin_popcountll(m); }
static inline unsigned MPOP(unsigned __int128 m) { return (unsigned)(__builtin_popcountll((unsigned long long)m) + __builtin_popcountll((unsigned long long)(m >> 64))); }
static inline unsigned MCTZ(uint32_t m) { return (unsigned)__builtin_ctz(m); }
static inline unsigned MCTZ(unsigned long long m) { return (unsigned)__builtin_ctzll(m); }
static inline unsigned MCTZ(unsigned __int128 m) { unsigned long long lo = (unsigned long long)m; return lo ? (unsigned)__builtin_ctzll(lo) : 64u + (unsigned)__builtin_ctzll((unsigned long long)(m >> 64)); }

void reverseBFS(uint32_t goal, vector<uint32_t>& dist) {
  dist.assign(N, NONE); dist[goal] = 0;
  vector<uint32_t> q; q.reserve(N); q.push_back(goal); size_t h = 0;
  while (h < q.size()) {
    uint32_t Lc = q[h++], d = dist[Lc];
    for (unsigned dir = 0; dir < 4; dir++) {
      uint32_t A = NB[Lc * 4 + OPP[dir]]; if (A == NONE || WALL[A]) continue;
      uint32_t K = NB[A * 4 + OPP[dir]]; if (K == NONE || WALL[K]) continue;
      if (dist[A] == NONE) { dist[A] = d + 1; q.push_back(A); }
    }
  }
}

void setup(const vector<string>& g) {
  R = g.size(); C = g[0].size(); N = R * C;
  WALL.assign(N, 0); ISGOAL.assign(N, 0); NB.assign(N * 4, NONE);
  goalCells.clear(); boxCells0.clear(); keeper0 = NONE;
  for (unsigned r = 0; r < R; r++) for (unsigned c = 0; c < C; c++) {
    unsigned v = (unsigned)(g[r][c] - '0'); uint32_t i = idx(r, c);
    if (v == 1) WALL[i] = 1;
    if (v == 4 || v == 5 || v == 6) { ISGOAL[i] = 1; goalCells.push_back(i); }
    if (v == 2 || v == 5) boxCells0.push_back(i);
    if (v == 3 || v == 6) keeper0 = i;
  }
  for (unsigned i = 0; i < N; i++) { int r = (int)(i / C), c = (int)(i % C); for (unsigned d = 0; d < 4; d++) { int nr = r + DR[d], nc = c + DC[d]; NB[i * 4 + d] = inb(nr, nc) ? idx((unsigned)nr, (unsigned)nc) : NONE; } }
  // dead squares via reverse-reachability from goals (NONE = max, so a plain
  // unsigned < is also the presence check)
  DMIN.assign(N, NONE);
  for (size_t gi = 0; gi < goalCells.size(); gi++) { vector<uint32_t> dist; reverseBFS(goalCells[gi], dist); for (unsigned i = 0; i < N; i++) if (dist[i] < DMIN[i]) DMIN[i] = dist[i]; }
  // live cells = non-dead (a box there can still reach a goal)
  liveIdx.assign(N, NONE); liveCell.clear();
  for (unsigned i = 0; i < N; i++) if (DMIN[i] != NONE) { liveIdx[i] = (uint32_t)liveCell.size(); liveCell.push_back(i); }
  L = (unsigned)liveCell.size();
}

inline bool deadCell(uint32_t c) { return DMIN[c] == NONE; }
template<class M> inline bool boxAt(M mask, uint32_t cell) { return cell != NONE && liveIdx[cell] != NONE && ((mask >> liveIdx[cell]) & (M)1); }

// heuristic: sum over boxes of push-distance to nearest goal (admissible for moves)
template<class M> inline uint32_t heuristic(M mask) {
  uint32_t s = 0; M m = mask;
  while (m) { unsigned b = MCTZ(m); m &= m - (M)1; s += DMIN[liveCell[b]]; }
  return s;
}

// freeze deadlock: box immovable on both axes (walls / dead / other frozen boxes)
template<class M> bool frozenRec(uint32_t cell, M mask, M assume);
template<class M> bool axisBlocked(uint32_t cell, unsigned axis, M mask, M assume) {
  uint32_t s1 = axis == 0 ? NB[cell * 4 + 2] : NB[cell * 4 + 0];
  uint32_t s2 = axis == 0 ? NB[cell * 4 + 3] : NB[cell * 4 + 1];
  if (s1 == NONE || WALL[s1] || s2 == NONE || WALL[s2]) return true;
  if (deadCell(s1) && deadCell(s2)) return true;
  uint32_t sides[2] = {s1, s2};
  for (unsigned k = 0; k < 2; k++) {
    uint32_t s = sides[k];
    if (liveIdx[s] != NONE && ((assume >> liveIdx[s]) & (M)1)) return true;
    if (boxAt(mask, s) && frozenRec<M>(s, mask, assume | ((M)1 << liveIdx[cell]))) return true;
  }
  return false;
}
template<class M> bool frozenRec(uint32_t cell, M mask, M assume) {
  return axisBlocked<M>(cell, 0, mask, assume) && axisBlocked<M>(cell, 1, mask, assume);
}
template<class M> inline bool isFrozen(uint32_t cell, M mask) { return frozenRec<M>(cell, mask, (M)0); }

// ---- keeper BFS distances (avoiding boxes). Per-thread scratch so the search
// can run on many cores at once ----
struct Scratch {
  vector<uint32_t> kdist, kstamp, queue; uint32_t kcur;
  uint64_t ex, gen;               // per-thread counters (no shared-cacheline atomics)
};
template<class M> void keeperBFS(Scratch& sc, M mask, uint32_t keeper) {
  if (sc.kcur == UINT32_MAX) { fill(sc.kstamp.begin(), sc.kstamp.end(), 0u); sc.kcur = 0; }  // stamp wrap guard
  sc.kcur++; sc.kstamp[keeper] = sc.kcur; sc.kdist[keeper] = 0;
  unsigned h = 0, tail = 0; sc.queue[tail++] = keeper;
  while (h < tail) {
    uint32_t cur = sc.queue[h++]; uint32_t base = cur * 4, d = sc.kdist[cur];
    for (unsigned dir = 0; dir < 4; dir++) {
      uint32_t n = NB[base + dir]; if (n == NONE || WALL[n] || boxAt(mask, n)) continue;
      if (sc.kstamp[n] != sc.kcur) { sc.kstamp[n] = sc.kcur; sc.kdist[n] = d + 1; sc.queue[tail++] = n; }
    }
  }
}
inline bool kReach(Scratch& sc, uint32_t cell) { return cell != NONE && sc.kstamp[cell] == sc.kcur; }

// ---- persistent worker pool: threads are spawned once and woken per wave via a
// generation counter (hundreds of waves would otherwise pay thread create+join
// each).  The mutex/cv handshake gives the happens-before for wave data both ways.
class Pool {
  unsigned nt = 0; vector<thread> ths; mutex mu; condition_variable cv, cvd;
  uint64_t gen = 0; unsigned running = 0; bool stop = false;
  function<void(unsigned)>* job = nullptr;
  void loop(unsigned wi) {
    uint64_t seen = 0;
    while (true) {
      function<void(unsigned)>* j;
      { unique_lock<mutex> lk(mu);
        cv.wait(lk, [&] { return stop || gen != seen; });
        if (stop) return;
        seen = gen; j = job; }
      (*j)(wi);
      { lock_guard<mutex> lk(mu); if (--running == 0) cvd.notify_all(); }
    }
  }
public:
  void init(unsigned n) { nt = n; for (unsigned i = 0; i < n; i++) ths.emplace_back([this, i] { loop(i); }); }
  void runAll(function<void(unsigned)>& j) {
    { lock_guard<mutex> lk(mu); job = &j; running = nt; gen++; }
    cv.notify_all();
    { unique_lock<mutex> lk(mu); cvd.wait(lk, [&] { return running == 0; }); }
  }
  void shutdown() {
    if (ths.empty()) return;
    { lock_guard<mutex> lk(mu); stop = true; }
    cv.notify_all();
    for (auto& t : ths) t.join();
    ths.clear();
  }
  ~Pool() { shutdown(); }
};

// ---- bucket entries: the fast path packs (mask,g,keeper) into one uint64
// (32+23+8 = 63 bits) to halve open-list memory traffic. The wide path uses a
// plain struct.  Uniform msk()/kpr()/gv() accessors keep the search generic. ----
template<class M> struct Ent {
  M mask; uint32_t keeper, g;
  static Ent make(M m, uint32_t k, uint32_t gg) { return Ent{m, k, gg}; }
  inline M msk() const { return mask; }
  inline uint32_t kpr() const { return keeper; }
  inline uint32_t gv() const { return g; }
};
template<> struct Ent<uint32_t> {
  uint64_t v;
  static Ent make(uint32_t m, uint32_t k, uint32_t gg) { return Ent{((uint64_t)m << 32) | ((uint64_t)gg << 8) | (uint8_t)k}; }
  inline uint32_t msk() const { return (uint32_t)(v >> 32); }
  inline uint32_t kpr() const { return (uint32_t)(v & 0xff); }
  inline uint32_t gv() const { return (uint32_t)((v >> 8) & GMASK); }
};

// ---- FAST PATH hash (<=32 live cells): lock-free open-addressed, ONE atomic word
// per slot: (key<<24) | value, key = (mask<<8)|keeper (<=40 bits).  All updates are
// atomic CAS, so many threads share it with no locks.  Fixed size (no rehash).
struct AHash {
  atomic<uint64_t>* a = nullptr; uint64_t capMask = 0;
  static constexpr uint64_t EMPTY = ~0ull;
  void init(size_t cap) {
    size_t s = 1; while (s < cap) s <<= 1;
    a = (atomic<uint64_t>*)bigalloc(s * sizeof(uint64_t)); capMask = s - 1;
    memset((void*)a, 0xFF, s * sizeof(uint64_t));          // all slots EMPTY
  }
  static inline uint64_t keyOf(uint32_t mask, uint32_t keeper) { return ((uint64_t)mask << 8) | keeper; }
  inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return x; }
  inline void prefetch(uint32_t mask, uint32_t keeper) { __builtin_prefetch(&a[mix(keyOf(mask, keeper)) & capMask]); }
  [[noreturn]] void full() { fprintf(stderr, "hash full -- increase hash size\n"); exit(2); }
  // insert or lower g. Returns true iff we set a new/better g on a non-closed slot.
  bool improve(uint32_t mask, uint32_t keeper, uint32_t ng) {
    uint64_t k = keyOf(mask, keeper); size_t i = mix(k) & capMask, probes = 0;
    while (true) {
      uint64_t cur = a[i].load(memory_order_relaxed);
      if (cur == EMPTY) { uint64_t des = (k << 24) | ng; if (a[i].compare_exchange_weak(cur, des, memory_order_relaxed)) return true; continue; }
      if ((cur >> 24) == k) {
        if (cur & CLOSEDBIT) return false;
        if ((uint32_t)(cur & GMASK) <= ng) return false;
        uint64_t des = (k << 24) | ng; if (a[i].compare_exchange_weak(cur, des, memory_order_relaxed)) return true; continue;
      }
      i = (i + 1) & capMask; if (++probes > capMask) full();
    }
  }
  // one-shot pop gate: state must exist with exactly g and not be closed. Close it.
  bool claim(uint32_t mask, uint32_t keeper, uint32_t g) {
    uint64_t k = keyOf(mask, keeper); size_t i = mix(k) & capMask, probes = 0;
    while (true) {
      uint64_t cur = a[i].load(memory_order_relaxed);
      if (cur == EMPTY) return false;
      if ((cur >> 24) == k) {
        uint64_t exp = (k << 24) | g;
        if (cur != exp) return false;                      // closed already, or a better g is known
        if (a[i].compare_exchange_strong(exp, exp | CLOSEDBIT, memory_order_relaxed)) return true;
        continue;                                          // raced. Re-inspect
      }
      i = (i + 1) & capMask; if (++probes > capMask) return false;
    }
  }
  size_t find(uint32_t mask, uint32_t keeper) { uint64_t k = keyOf(mask, keeper); size_t i = mix(k) & capMask, probes = 0; while (true) { uint64_t cur = a[i].load(memory_order_relaxed); if (cur == EMPTY) return NPOS; if ((cur >> 24) == k) return i; i = (i + 1) & capMask; if (++probes > capMask) return NPOS; } }
  inline uint32_t val(size_t i) { return (uint32_t)(a[i].load(memory_order_relaxed) & VALMASK); }
  size_t slots() { return capMask + 1; }
  bool occupied(size_t i, uint32_t& mask, uint32_t& keeper, uint32_t& v) { uint64_t cur = a[i].load(memory_order_relaxed); if (cur == EMPTY) return false; uint64_t key = cur >> 24; mask = (uint32_t)(key >> 8); keeper = (uint32_t)(key & 0xff); v = (uint32_t)(cur & VALMASK); return true; }
};

// ---- WIDE PATH hash (33..128 live cells): the mask no longer fits beside the
// value in one word, so each slot is an atomic control word (EMPTY / WRITING /
// READY|value) plus side arrays holding the full key.  A slot is claimed EMPTY->
// WRITING, its key written, then published READY with release/acquire ordering so
// readers that observe READY also observe the key.  Readers spin briefly on a
// slot mid-publication, so this is not formally lock-free (unlike AHash).
template<class M>
struct WideHash {
  atomic<uint64_t>* ctrl = nullptr; M* km = nullptr; uint8_t* kk = nullptr; uint64_t capMask = 0;
  static constexpr uint64_t EMPTY = 0ull;
  static constexpr uint64_t WRITING = 1ull << 62;
  static constexpr uint64_t READY = 1ull << 63;
  void init(size_t cap) {
    size_t s = 1; while (s < cap) s <<= 1;
    ctrl = (atomic<uint64_t>*)bigalloc(s * sizeof(uint64_t));
    km = (M*)bigalloc(s * sizeof(M)); kk = (uint8_t*)bigalloc(s); capMask = s - 1;
    memset((void*)ctrl, 0, s * sizeof(uint64_t));          // all slots EMPTY
  }
  inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return x; }
  inline uint64_t hashKey(M mask, uint32_t keeper) {
    uint64_t h = (uint64_t)mask;
    if constexpr (sizeof(M) > 8) h ^= mix((uint64_t)(mask >> 64));
    return mix(h ^ ((uint64_t)keeper * 0x9E3779B97F4A7C15ull));
  }
  inline void prefetch(M mask, uint32_t keeper) { __builtin_prefetch(&ctrl[hashKey(mask, keeper) & capMask]); }
  [[noreturn]] void full() { fprintf(stderr, "WideHash full -- increase hash size\n"); exit(2); }
  bool improve(M mask, uint32_t keeper, uint32_t ng) {
    size_t i = hashKey(mask, keeper) & capMask, probes = 0;
    while (true) {
      uint64_t c = ctrl[i].load(memory_order_acquire);
      if (c == EMPTY) {
        uint64_t exp = EMPTY;
        if (ctrl[i].compare_exchange_strong(exp, WRITING, memory_order_acq_rel)) {
          km[i] = mask; kk[i] = (uint8_t)keeper;
          ctrl[i].store(READY | ng, memory_order_release);
          return true;
        }
        continue;                                              // lost the claim. Re-inspect this slot
      }
      if (c == WRITING) continue;                              // another thread is publishing. Spin
      if (km[i] == mask && kk[i] == (uint8_t)keeper) {         // key fixed once READY -> safe to read
        uint32_t v = (uint32_t)(c & VALMASK);
        if (v & CLOSEDBIT) return false;
        if ((v & GMASK) <= ng) return false;
        if (ctrl[i].compare_exchange_strong(c, READY | ng, memory_order_acq_rel)) return true;
        continue;
      }
      i = (i + 1) & capMask; if (++probes > capMask) full();
    }
  }
  // one-shot pop gate: state must exist with exactly g and not be closed. Close it.
  bool claim(M mask, uint32_t keeper, uint32_t g) {
    size_t i = hashKey(mask, keeper) & capMask, probes = 0;
    while (true) {
      uint64_t c = ctrl[i].load(memory_order_acquire);
      if (c == EMPTY) return false;
      if (c == WRITING) continue;
      if (km[i] == mask && kk[i] == (uint8_t)keeper) {
        uint64_t exp = READY | g;
        if (c != exp) return false;                            // closed already, or a better g is known
        if (ctrl[i].compare_exchange_strong(exp, exp | CLOSEDBIT, memory_order_acq_rel)) return true;
        continue;                                              // raced. Re-inspect
      }
      i = (i + 1) & capMask; if (++probes > capMask) return false;
    }
  }
  size_t find(M mask, uint32_t keeper) {
    size_t i = hashKey(mask, keeper) & capMask, probes = 0;
    while (true) {
      uint64_t c = ctrl[i].load(memory_order_acquire);
      if (c == EMPTY) return NPOS;
      if (c == WRITING) continue;
      if (km[i] == mask && kk[i] == (uint8_t)keeper) return i;
      i = (i + 1) & capMask; if (++probes > capMask) return NPOS;
    }
  }
  inline uint32_t val(size_t i) { return (uint32_t)(ctrl[i].load(memory_order_acquire) & VALMASK); }
  size_t slots() { return capMask + 1; }
  bool occupied(size_t i, M& mask, uint32_t& keeper, uint32_t& v) { uint64_t c = ctrl[i].load(memory_order_acquire); if (!(c & READY)) return false; mask = km[i]; keeper = kk[i]; v = (uint32_t)(c & VALMASK); return true; }
};

// ---- helpers for solution reconstruction (single-threaded, run once at the end) ----
// keeper BFS distances in a given box mask (NONE = unreachable). Avoids walls/boxes.
template<class M> static vector<uint32_t> keeperDist(M mask, uint32_t src) {
  vector<uint32_t> dist(N, NONE); if (src == NONE) return dist; dist[src] = 0;
  vector<uint32_t> q; q.reserve(N); q.push_back(src); size_t h = 0;
  while (h < q.size()) {
    uint32_t cur = q[h++], base = cur * 4, dc = dist[cur];
    for (unsigned d = 0; d < 4; d++) { uint32_t n = NB[base + d]; if (n == NONE || WALL[n] || boxAt(mask, n)) continue; if (dist[n] == NONE) { dist[n] = dc + 1; q.push_back(n); } }
  }
  return dist;
}
// shortest keeper walk from src to dst as a lowercase move string ("udlr"), "" if src==dst.
template<class M> static string keeperPathMoves(M mask, uint32_t src, uint32_t dst) {
  if (src == dst) return string();
  vector<uint32_t> par(N, NONE), pdir(N, NONE); vector<char> vis(N, 0);
  vector<uint32_t> q; q.push_back(src); vis[src] = 1; size_t h = 0; bool done = false;
  while (h < q.size() && !done) {
    uint32_t cur = q[h++], base = cur * 4;
    for (unsigned d = 0; d < 4; d++) { uint32_t n = NB[base + d]; if (n == NONE || WALL[n] || boxAt(mask, n) || vis[n]) continue; vis[n] = 1; par[n] = cur; pdir[n] = d; if (n == dst) { done = true; break; } q.push_back(n); }
  }
  string s; uint32_t cur = dst; while (cur != src && par[cur] != NONE) { s += "udlr"[pdir[cur]]; cur = par[cur]; }
  reverse(s.begin(), s.end()); return s;
}

// ---- the search + optional reconstruction, templated on the box-mask type M and
// its matching hash H (AHash for uint32, WideHash<M> for the wide case) ----
template<class M, class H>
int run(const string& board, double timeout, bool wantSolution, bool dbg) {
  // a start box on a dead square (liveIdx == NONE) can never reach any goal.
  // catch it BEFORE building the mask. A NONE shift would be UB and the box
  // would vanish.
  for (uint32_t bc : boxCells0) if (liveIdx[bc] == NONE) { printf("RESULT %s optimal=UNSOLVABLE\n", board.c_str()); return 0; }
  M boxMask0 = 0; for (uint32_t bc : boxCells0) boxMask0 |= (M)1 << liveIdx[bc];
  M goalBoxMask = 0; for (uint32_t gc : goalCells) if (liveIdx[gc] != NONE) goalBoxMask |= (M)1 << liveIdx[gc];

  unsigned nb = MPOP(boxMask0);
  printf("[%s] boxes=%u goals=%zu live=%u h0=%u timeout=%.0fs\n", board.c_str(), nb, goalCells.size(), L, heuristic(boxMask0), timeout);
  fflush(stdout);

  if (dbg) {
    M m = boxMask0; uint32_t tot = 0; while (m) { unsigned b = MCTZ(m); m &= m - (M)1; uint32_t cell = liveCell[b]; printf("box cell %u (r%u,c%u) DMIN=%u\n", cell, cell / C, cell % C, DMIN[cell]); tot += DMIN[cell]; } printf("C++ SUM h0=%u\n", tot);
    uint32_t gcell = 6 * C + 6;
    if (gcell < N) {
      vector<uint32_t> dd; reverseBFS(gcell, dd);
      printf("reverseBFS(66=%u): dist[65]=%d dist[56]=%d dist[66]=%d\n", gcell, (int)dd[6 * C + 5], (int)dd[5 * C + 6], (int)dd[gcell]);
      printf("ISGOAL[66]=%d NB[66] U=%d D=%d L=%d R=%d\n", (int)ISGOAL[gcell], (int)NB[gcell * 4 + 0], (int)NB[gcell * 4 + 1], (int)NB[gcell * 4 + 2], (int)NB[gcell * 4 + 3]);
    }
    return 0;
  }

  // ---- parallel move-optimal A* (level-synchronous Dial) ----
  unsigned NT = thread::hardware_concurrency(); if (NT == 0) NT = 4; if (NT > 32) NT = 32;
  // hash sizing: fast path keeps the original fixed sizes. Wide path sizes from a
  // reachable-state estimate (few boxes on a big board) capped for memory.
  size_t hashSize;
  if constexpr (sizeof(M) <= 4) {
    hashSize = (L >= 24) ? ((size_t)1 << 27) : ((size_t)1 << 22);
  } else {
    long double comb = 1;
    for (unsigned i = 0; i < nb; i++) { comb *= (long double)(L - i) / (long double)(i + 1); if (comb > 1e13L) { comb = 1e13L; break; } }
    long double est = comb * (long double)N * 3.0L;
    long double cap = (long double)((size_t)1 << 25);
    size_t want = (size_t)(est < cap ? est : cap);
    hashSize = 1 << 16; while (hashSize < want) hashSize <<= 1;
  }
  H hash; hash.init(hashSize);

  constexpr uint32_t MAXF = 6000;
  using E = Ent<M>;
  vector<vector<vector<E>>> bucket(MAXF);           // per f: list of chunks (bulk-moved, never copied)

  uint32_t h0 = heuristic(boxMask0);
  if (h0 >= MAXF) { printf("RESULT %s optimal=UNSUPPORTED (h0=%u exceeds MAXF=%u)\n", board.c_str(), h0, MAXF); return 0; }
  hash.improve(boxMask0, keeper0, 0);
  bucket[h0].push_back(vector<E>{E::make(boxMask0, keeper0, 0)});

  atomic<uint32_t> best(NONE);                      // NONE = no solution known yet
  Pool pool; bool poolInited = false;               // spawned on first parallel wave
  vector<Scratch> scratch(NT);
  for (auto& sc : scratch) { sc.kdist.assign(N, 0); sc.kstamp.assign(N, 0); sc.queue.assign(N, 0); sc.kcur = 0; sc.ex = sc.gen = 0; }
  // per-(thread,f) successor chunks + touched-f lists: workers append locally,
  // the merge step moves whole chunks. No per-element copying, no locks.
  vector<vector<vector<E>>> tf(NT, vector<vector<E>>(MAXF));
  vector<vector<uint32_t>> touched(NT);
  struct Cand { M nmask; uint32_t bcell, ng, nf; };
  vector<vector<Cand>> candBuf(NT); for (auto& cb : candBuf) cb.reserve(256);
  uint64_t lastReport = 0;
  bool timedOut = false; uint32_t lbAtTimeout = 0;
  auto t0 = chrono::steady_clock::now();
  auto sumEx = [&]() { uint64_t s = 0; for (auto& sc : scratch) s += sc.ex; return s; };
  auto sumGen = [&]() { uint64_t s = 0; for (auto& sc : scratch) s += sc.gen; return s; };

  // expand one wave of same-f entries. Returns the successors that land back in f.
  auto processWave = [&](vector<E>& wave, uint32_t f) -> vector<E> {
    auto worker = [&](unsigned wi, size_t stride) {
      Scratch& sc = scratch[wi];
      auto& myTf = tf[wi]; auto& myTouched = touched[wi]; auto& cands = candBuf[wi];
      uint64_t ex = 0, gen = 0;
      for (size_t j = wi; j < wave.size(); j += stride) {
        if (j + stride < wave.size()) { const E& e2 = wave[j + stride]; hash.prefetch(e2.msk(), e2.kpr()); }
        M mask = wave[j].msk(); uint32_t keeper = wave[j].kpr(), g = wave[j].gv();
        if (!hash.claim(mask, keeper, g)) continue;            // stale, closed, or another thread won it
        ex++;
        keeperBFS(sc, mask, keeper);
        // all boxes placed? -> finish by walking keeper to nearest free goal
        // (NONE = max, so the min-update needs no separate presence check)
        if ((mask & goalBoxMask) == mask) {
          uint32_t fw = NONE;
          for (uint32_t gc : goalCells) if (!boxAt(mask, gc) && kReach(sc, gc) && sc.kdist[gc] < fw) fw = sc.kdist[gc];
          if (fw != NONE) { uint32_t cand = g + fw; uint32_t cb = best.load(memory_order_relaxed); while (cand < cb && !best.compare_exchange_weak(cb, cand, memory_order_relaxed)); }
        }
        // pass 1: collect legal pushes, prefetch their hash slots.  h(mask) = f - g
        // exactly (entries are bucketed by f = g + h), so successor f is an O(1)
        // delta instead of a loop over boxes.
        uint32_t hpar = f - g;
        cands.clear();
        M m = mask;
        while (m) {
          unsigned b = MCTZ(m); m &= m - (M)1; uint32_t bcell = liveCell[b];
          for (unsigned d = 0; d < 4; d++) {
            uint32_t pf = NB[bcell * 4 + OPP[d]]; if (!kReach(sc, pf)) continue;
            uint32_t t = NB[bcell * 4 + d];
            if (t == NONE || WALL[t] || deadCell(t) || boxAt(mask, t)) continue;
            M nmask = (mask & ~((M)1 << b)) | ((M)1 << liveIdx[t]);
            if (!ISGOAL[t] && isFrozen(t, nmask)) continue;
            uint32_t ng = g + sc.kdist[pf] + 1;
            uint32_t nf = ng + hpar - DMIN[bcell] + DMIN[t];
#ifdef HCHECK
            if (nf != ng + heuristic(nmask)) { fprintf(stderr, "HCHECK FAIL f=%u g=%u\n", f, g); abort(); }
#endif
            hash.prefetch(nmask, bcell);
            cands.push_back({nmask, bcell, ng, nf});
          }
        }
        // pass 2: relax into the hash (slots are now in flight).  Successors with
        // nf >= best can never be expanded (levels stop at best, and best only
        // decreases), so skip queueing them. But the improve() must still run so
        // the g-map stays complete for solution reconstruction.
        uint32_t cb = best.load(memory_order_relaxed);
        for (const Cand& cd : cands) {
          if (hash.improve(cd.nmask, cd.bcell, cd.ng)) {
            if (cd.nf < MAXF && cd.nf < cb) {
              auto& dst = myTf[cd.nf];
              if (dst.empty()) myTouched.push_back(cd.nf);
              dst.push_back(E::make(cd.nmask, cd.bcell, cd.ng));
              gen++;
            }
          }
        }
      }
      sc.ex += ex; sc.gen += gen;
    };
    if (NT == 1 || wave.size() < 2000) { worker(0, 1); }           // small wave: skip pool overhead
    else {
      if (!poolInited) { pool.init(NT); poolInited = true; }
      function<void(unsigned)> j = [&](unsigned wi) { worker(wi, (size_t)NT); };
      pool.runAll(j);
    }
    // merge: same-f chunks become the next wave. Cross-f chunks are MOVED into
    // their bucket (single-threaded, pointer moves only. No element copies).
    size_t nextTotal = 0;
    for (unsigned t = 0; t < NT; t++) for (uint32_t ff : touched[t]) if (ff == f) nextTotal += tf[t][ff].size();
    vector<E> next; next.reserve(nextTotal);
    for (unsigned t = 0; t < NT; t++) {
      for (uint32_t ff : touched[t]) {
        auto& ch = tf[t][ff];
        if (ff == f) { next.insert(next.end(), ch.begin(), ch.end()); ch.clear(); }
        else { bucket[ff].push_back(move(ch)); ch = vector<E>(); }
      }
      touched[t].clear();
    }
    return next;
  };

  for (uint32_t f = 0; f < MAXF && !timedOut; f++) {
    if (f >= best.load(memory_order_relaxed)) break;   // admissible: nothing below best remains
    if (bucket[f].empty()) continue;
    size_t tot = 0; for (auto& ch : bucket[f]) tot += ch.size();
    vector<E> wave; wave.reserve(tot);
    for (auto& ch : bucket[f]) wave.insert(wave.end(), ch.begin(), ch.end());
    bucket[f].clear(); bucket[f].shrink_to_fit();
    while (!wave.empty()) {
      wave = processWave(wave, f);
      // once best <= f, every cheaper level is exhausted, so best is provably
      // optimal. The rest of this level can only rediscover it.
      if (best.load(memory_order_relaxed) <= f) { wave.clear(); break; }
      uint64_t ex = sumEx();
      if (ex - lastReport >= 4000000) {
        lastReport = ex;
        double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        uint32_t bb = best.load(memory_order_relaxed);
        printf("  ...expanded %llu, f=%u gen=%llu best=%d %.0fs\n", (unsigned long long)ex, f, (unsigned long long)sumGen(), bb == NONE ? -1 : (int)bb, el);
        fflush(stdout);
        if (el > timeout) { timedOut = true; lbAtTimeout = f; break; }
      }
    }
  }

  double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
  uint32_t bestv = best.load(memory_order_relaxed);
  uint64_t ex = sumEx(), gen = sumGen();
  if (timedOut) { printf("RESULT %s optimal=TIMEOUT lower_bound=%u best_found=%s expanded=%llu seconds=%.1f threads=%u\n", board.c_str(), lbAtTimeout, bestv == NONE ? "none" : to_string(bestv).c_str(), (unsigned long long)ex, el, NT); return 0; }
  if (bestv != NONE) printf("RESULT %s optimal=%u expanded=%llu generated=%llu seconds=%.1f threads=%u\n", board.c_str(), bestv, (unsigned long long)ex, (unsigned long long)gen, el, NT);
  else printf("RESULT %s optimal=NONE expanded=%llu seconds=%.1f threads=%u\n", board.c_str(), (unsigned long long)ex, el, NT);

  // ---- optional: reconstruct one optimal move path from the completed g-map ----
  // The hash still holds the optimal g of every explored (boxMask, keeper) state
  // (closing only sets a high bit. g survives in the low bits).  We (1) find a
  // finishing state F. All boxes on goals with g(F)+walk(keeper->free goal)==opt
  // . Then (2) walk backwards: at each state the keeper sits on the cell the last
  // pushed box vacated, so its just-moved box is an orthogonal neighbour. We invert
  // that push to a predecessor whose stored g satisfies g_pred + (keeper walk) + 1
  // == g.  Every explored state's g is its true distance from start, so any g-tight
  // predecessor provably chains back to the start. Greedy first-match needs no
  // backtracking.  Reversing the collected macro-moves and appending the finishing
  // walk yields a full optimal path.
  if (wantSolution && bestv != NONE) {
    // (1) locate a finishing state F consistent with the optimum.
    M maskF = 0; uint32_t keeperF = NONE, freeGoalF = NONE;
    for (size_t i = 0; i < hash.slots(); i++) {
      M mask; uint32_t keeper, v;
      if (!hash.occupied(i, mask, keeper, v)) continue;
      uint32_t g = v & GMASK;
      if ((mask & goalBoxMask) != mask) continue;                    // not all boxes on goals
      vector<uint32_t> dist = keeperDist(mask, keeper);
      uint32_t fw = NONE, fg = NONE;
      for (uint32_t gc : goalCells) if (!boxAt(mask, gc) && dist[gc] < fw) { fw = dist[gc]; fg = gc; }
      if (fw != NONE && g + fw == bestv) { maskF = mask; keeperF = keeper; freeGoalF = fg; break; }
    }
    if (keeperF == NONE) { printf("SOLUTION %s reconstruction_failed (no finishing state)\n", board.c_str()); return 0; }

    string finish = keeperPathMoves(maskF, keeperF, freeGoalF);       // tail walk onto the free goal
    // (2) invert pushes back to the start state.
    vector<string> segsRev; M mask = maskF; uint32_t keeper = keeperF; bool ok = true;
    while (!(mask == boxMask0 && keeper == keeper0)) {
      size_t vi = hash.find(mask, keeper); if (vi == NPOS) { ok = false; break; }
      uint32_t g = hash.val(vi) & GMASK; bool found = false;
      for (unsigned d = 0; d < 4 && !found; d++) {
        uint32_t t = NB[keeper * 4 + d]; if (!boxAt(mask, t)) continue;          // box that was just pushed
        uint32_t pf = NB[keeper * 4 + OPP[d]]; if (pf == NONE || WALL[pf]) continue;  // keeper stood here to push
        M pmask = (mask & ~((M)1 << liveIdx[t])) | ((M)1 << liveIdx[keeper]);    // undo push: box back to keeper cell
        if (boxAt(pmask, pf)) continue;                                          // push spot must have been free
        vector<uint32_t> dist = keeperDist(pmask, pf);
        for (uint32_t kp = 0; kp < N; kp++) {
          if (dist[kp] == NONE) continue;
          size_t vj = hash.find(pmask, kp); if (vj == NPOS) continue;
          uint32_t gp = hash.val(vj) & GMASK;
          if (gp + dist[kp] + 1 == g) {                                          // g-tight predecessor
            string seg = keeperPathMoves(pmask, kp, pf); seg += "UDLR"[d];
            segsRev.push_back(seg); mask = pmask; keeper = kp; found = true; break;
          }
        }
      }
      if (!found) { ok = false; break; }
    }
    if (!ok) { printf("SOLUTION %s reconstruction_failed\n", board.c_str()); return 0; }
    string sol; for (size_t i = segsRev.size(); i-- > 0; ) sol += segsRev[i]; sol += finish;
    printf("SOLUTION %s moves=%zu %s\n", board.c_str(), sol.size(), sol.c_str());
    if ((uint32_t)sol.size() != bestv) printf("WARNING %s solution length %zu != optimal %u\n", board.c_str(), sol.size(), bestv);
    fflush(stdout);
  }
  return 0;
}

int main(int argc, char** argv) {          // argc/argv/return are int by the language
  // Parse args: positionals are <boardId> [timeoutSec]. Flags may appear anywhere.
  string board = "p15"; double timeout = 600.0; bool wantSolution = false, dbg = false;
  {
    vector<string> pos;
    for (int i = 1; i < argc; i++) {
      string a = argv[i];
      if (a == "--solution" || a == "-s" || a == "sol") wantSolution = true;
      else if (a == "dbg") dbg = true;
      else pos.push_back(a);
    }
    if (pos.size() >= 1) board = pos[0];
    if (pos.size() >= 2) timeout = atof(pos[1].c_str());
  }
  if (!BOARDS.count(board)) { printf("unknown board %s\n", board.c_str()); return 1; }
  setup(BOARDS[board]);

  // Dispatch on live-cell count: <=32 uses the fast uint32 mask + packed hash.
  // 33..128 uses a 128-bit mask + wide hash. Beyond that is unsupported.
  if (N > 255) { printf("RESULT %s optimal=UNSUPPORTED (%u cells exceeds 8-bit keeper index)\n", board.c_str(), N); return 0; }
  if (L <= 32) return run<uint32_t, AHash>(board, timeout, wantSolution, dbg);
  if (L <= 128) return run<unsigned __int128, WideHash<unsigned __int128>>(board, timeout, wantSolution, dbg);
  printf("RESULT %s optimal=UNSUPPORTED (live=%u exceeds the 128-bit box mask)\n", board.c_str(), L);
  return 0;
}
