// Move-OPTIMAL Sokoban solver (CS161 variant: keeper must also end on a goal).
//
// Techniques:
//  * Compact encoding: boxes = a bitmask over the <=32 "live" cells (uint32),
//    keeper = a cell index (uint8).  A whole state key is one uint64
//    ((boxMask<<8)|keeper) -> ~16 bytes in the hash vs Python's ~1.3 KB.
//  * Open-addressed flat hash (linear probing) for best-g / closed.
//  * Dial bucket queue for the A* open list (integer f) -- entries are just the
//    64-bit state key; g is recovered as f - h(mask) since the heuristic depends
//    only on the box mask.
//  * Macro-push moves with move-cost edges: each edge is (shortest keeper walk
//    to the push spot) + 1, so g counts real keeper MOVES and the answer is the
//    true minimum-move solution, while the state space stays (boxes, keeper).
//  * Dead-square table (reverse-reachability from goals) + freeze-deadlock
//    pruning to cut hopeless branches.
//
// Reports the optimal move count (not the path).  Usage: ./optimal <boardId> [timeoutSec]

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <chrono>
#include <thread>
#include <atomic>
#include <algorithm>
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

int R, C, N;
vector<uint8_t> WALL, ISGOAL;
vector<int> NB;                 // N*4 neighbor cell or -1
vector<int> DMIN;               // min push-distance to any goal, -1 = dead
vector<int> liveIdx;            // cell -> bit index (0..L-1) or -1
vector<int> liveCell;           // bit index -> cell
vector<int> goalCells;
int L;                          // number of live cells
uint32_t boxMask0; int keeper0;
const int DR[4] = {-1, 1, 0, 0}, DC[4] = {0, 0, -1, 1}, OPP[4] = {1, 0, 3, 2};

inline int idx(int r, int c) { return r * C + c; }
inline bool inb(int r, int c) { return r >= 0 && r < R && c >= 0 && c < C; }

void reverseBFS(int goal, vector<int>& dist) {
  dist.assign(N, -1); dist[goal] = 0;
  vector<int> q; q.reserve(N); q.push_back(goal); size_t h = 0;
  while (h < q.size()) {
    int Lc = q[h++], d = dist[Lc];
    for (int dir = 0; dir < 4; dir++) {
      int A = NB[Lc * 4 + OPP[dir]]; if (A < 0 || WALL[A]) continue;
      int K = NB[A * 4 + OPP[dir]]; if (K < 0 || WALL[K]) continue;
      if (dist[A] < 0) { dist[A] = d + 1; q.push_back(A); }
    }
  }
}

void setup(const vector<string>& g) {
  R = g.size(); C = g[0].size(); N = R * C;
  WALL.assign(N, 0); ISGOAL.assign(N, 0); NB.assign(N * 4, -1);
  goalCells.clear(); boxMask0 = 0; keeper0 = -1;
  vector<int> boxCells;
  for (int r = 0; r < R; r++) for (int c = 0; c < C; c++) {
    int v = g[r][c] - '0', i = idx(r, c);
    if (v == 1) WALL[i] = 1;
    if (v == 4 || v == 5 || v == 6) { ISGOAL[i] = 1; goalCells.push_back(i); }
    if (v == 2 || v == 5) boxCells.push_back(i);
    if (v == 3 || v == 6) keeper0 = i;
  }
  for (int i = 0; i < N; i++) { int r = i / C, c = i % C; for (int d = 0; d < 4; d++) { int nr = r + DR[d], nc = c + DC[d]; NB[i * 4 + d] = inb(nr, nc) ? idx(nr, nc) : -1; } }
  // dead squares via reverse-reachability from goals
  DMIN.assign(N, -1);
  for (int gi = 0; gi < (int)goalCells.size(); gi++) { vector<int> dist; reverseBFS(goalCells[gi], dist); for (int i = 0; i < N; i++) if (dist[i] >= 0 && (DMIN[i] < 0 || dist[i] < DMIN[i])) DMIN[i] = dist[i]; }
  // live cells = non-dead (a box there can still reach a goal)
  liveIdx.assign(N, -1); liveCell.clear();
  for (int i = 0; i < N; i++) if (DMIN[i] >= 0) { liveIdx[i] = liveCell.size(); liveCell.push_back(i); }
  L = liveCell.size();
  for (int bc : boxCells) boxMask0 |= (1u << liveIdx[bc]);
}

inline bool deadCell(int c) { return DMIN[c] < 0; }
inline bool boxAt(uint32_t mask, int cell) { return cell >= 0 && liveIdx[cell] >= 0 && ((mask >> liveIdx[cell]) & 1u); }

// heuristic: sum over boxes of push-distance to nearest goal (admissible for moves)
inline int heuristic(uint32_t mask) {
  int s = 0; uint32_t m = mask;
  while (m) { int b = __builtin_ctz(m); m &= m - 1; s += DMIN[liveCell[b]]; }
  return s;
}

// freeze deadlock: box immovable on both axes (walls / dead / other frozen boxes)
bool frozenRec(int cell, uint32_t mask, uint64_t assume);
bool axisBlocked(int cell, int axis, uint32_t mask, uint64_t assume) {
  int s1 = axis == 0 ? NB[cell * 4 + 2] : NB[cell * 4 + 0];
  int s2 = axis == 0 ? NB[cell * 4 + 3] : NB[cell * 4 + 1];
  if (s1 < 0 || WALL[s1] || s2 < 0 || WALL[s2]) return true;
  if (deadCell(s1) && deadCell(s2)) return true;
  int sides[2] = {s1, s2};
  for (int k = 0; k < 2; k++) {
    int s = sides[k];
    if (liveIdx[s] >= 0 && (assume >> liveIdx[s]) & 1ull) return true;
    if (boxAt(mask, s) && frozenRec(s, mask, assume | (1ull << liveIdx[cell]))) return true;
  }
  return false;
}
bool frozenRec(int cell, uint32_t mask, uint64_t assume) {
  return axisBlocked(cell, 0, mask, assume) && axisBlocked(cell, 1, mask, assume);
}
inline bool isFrozen(int cell, uint32_t mask) { return frozenRec(cell, mask, 0); }

// ---- keeper BFS distances (avoiding boxes); per-thread scratch so the search
// can run on many cores at once ----
struct Scratch { vector<int> kdist, kstamp, queue; int kcur; };
void keeperBFS(Scratch& sc, uint32_t mask, int keeper) {
  sc.kcur++; sc.kstamp[keeper] = sc.kcur; sc.kdist[keeper] = 0;
  int h = 0, tail = 0; sc.queue[tail++] = keeper;
  while (h < tail) {
    int cur = sc.queue[h++]; int base = cur * 4, d = sc.kdist[cur];
    for (int dir = 0; dir < 4; dir++) {
      int n = NB[base + dir]; if (n < 0 || WALL[n] || boxAt(mask, n)) continue;
      if (sc.kstamp[n] != sc.kcur) { sc.kstamp[n] = sc.kcur; sc.kdist[n] = d + 1; sc.queue[tail++] = n; }
    }
  }
}
inline bool kReach(Scratch& sc, int cell) { return cell >= 0 && sc.kstamp[cell] == sc.kcur; }

// ---- lock-free open-addressed hash, ONE atomic word per slot:
// (key<<24) | value, value = (closed at bit 23) | (best g in bits 0..22).
// The key is <=40 bits so it never overlaps the value, EMPTY = ~0 never occurs
// as a real slot.  All updates are atomic CAS, so many threads share it with no
// locks.  Fixed size (no rehash while parallel) -- sized up front for the board.
struct AHash {
  atomic<uint64_t>* a = nullptr; uint64_t capMask = 0;
  static constexpr uint64_t EMPTY = ~0ull;
  static constexpr uint32_t VALMASK = (1u << 24) - 1;
  static constexpr uint32_t CLOSED = 1u << 23;
  void init(size_t cap) {
    size_t s = 1; while (s < cap) s <<= 1;
    a = new atomic<uint64_t>[s]; capMask = s - 1;
    for (size_t i = 0; i < s; i++) a[i].store(EMPTY, memory_order_relaxed);
  }
  inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return x; }
  // insert or lower g; returns true iff we set a new/better g on a non-closed slot.
  bool improve(uint64_t k, uint32_t ng) {
    size_t i = mix(k) & capMask;
    while (true) {
      uint64_t cur = a[i].load(memory_order_relaxed);
      if (cur == EMPTY) { uint64_t des = (k << 24) | ng; if (a[i].compare_exchange_weak(cur, des, memory_order_relaxed)) return true; continue; }
      if ((cur >> 24) == k) {
        if (cur & CLOSED) return false;
        if ((uint32_t)(cur & 0x7fffff) <= ng) return false;
        uint64_t des = (k << 24) | ng; if (a[i].compare_exchange_weak(cur, des, memory_order_relaxed)) return true; continue;
      }
      i = (i + 1) & capMask;
    }
  }
  long find(uint64_t k) { size_t i = mix(k) & capMask; while (true) { uint64_t cur = a[i].load(memory_order_relaxed); if (cur == EMPTY) return -1; if ((cur >> 24) == k) return (long)i; i = (i + 1) & capMask; } }
  inline uint32_t val(long i) { return (uint32_t)(a[i].load(memory_order_relaxed) & VALMASK); }
  // try to atomically mark slot (key k, current g==expg, not closed) as closed.
  bool tryClose(long i, uint64_t k, uint32_t expg) { uint64_t exp = (k << 24) | expg; return a[i].compare_exchange_strong(exp, exp | CLOSED, memory_order_relaxed); }
};

int main(int argc, char** argv) {
  string board = argc > 1 ? argv[1] : "p15";
  double timeout = argc > 2 ? atof(argv[2]) : 600.0;
  if (!BOARDS.count(board)) { printf("unknown board %s\n", board.c_str()); return 1; }
  setup(BOARDS[board]);
  uint32_t goalBoxMask = 0; for (int gc : goalCells) if (liveIdx[gc] >= 0) goalBoxMask |= (1u << liveIdx[gc]);

  int nb = __builtin_popcount(boxMask0);
  printf("[%s] boxes=%d goals=%d live=%d h0=%d timeout=%.0fs\n", board.c_str(), nb, (int)goalCells.size(), L, heuristic(boxMask0), timeout);
  fflush(stdout);
  if (L > 32) { printf("RESULT %s optimal=UNSUPPORTED (live=%d exceeds the 32-bit box mask)\n", board.c_str(), L); return 0; }
  if (argc > 3 && string(argv[3]) == "dbg") {
    uint32_t m = boxMask0; int tot = 0; while (m) { int b = __builtin_ctz(m); m &= m - 1; int cell = liveCell[b]; printf("box cell %d (r%d,c%d) DMIN=%d\n", cell, cell / C, cell % C, DMIN[cell]); tot += DMIN[cell]; } printf("C++ SUM h0=%d\n", tot);
    int gcell = 6 * C + 6;
    vector<int> dd; reverseBFS(gcell, dd);
    printf("reverseBFS(66=%d): dist[65]=%d dist[56]=%d dist[66]=%d\n", gcell, dd[6*C+5], dd[5*C+6], dd[gcell]);
    printf("ISGOAL[66]=%d NB[66] U=%d D=%d L=%d R=%d\n", ISGOAL[gcell], NB[gcell*4+0], NB[gcell*4+1], NB[gcell*4+2], NB[gcell*4+3]);
    printf("WALL[65]=%d WALL[64]=%d WALL[56]=%d WALL[46]=%d\n", (int)WALL[6*C+5], (int)WALL[6*C+4], (int)WALL[5*C+6], (int)WALL[4*C+6]);
    return 0;
  }

  // any box already on a dead square? unsolvable
  {
    uint32_t m = boxMask0; while (m) { int b = __builtin_ctz(m); m &= m - 1; if (deadCell(liveCell[b]) && !ISGOAL[liveCell[b]]) { printf("RESULT %s optimal=UNSOLVABLE\n", board.c_str()); return 0; } }
  }

  // ---- parallel move-optimal A* (level-synchronous Dial) ----
  // The open list is bucketed by f.  We process levels in increasing f; within a
  // level, entries are expanded in "waves" by a pool of worker threads.  Every
  // shared write (best-g map, closed bit, global best) goes through an atomic, so
  // there are no locks.  Because the heuristic is consistent, the first time a
  // state is popped it is popped at its minimum g, so closing on first pop stays
  // optimal even across threads.
  unsigned NT = thread::hardware_concurrency(); if (NT == 0) NT = 4; if (NT > 14) NT = 14;
  size_t hashSize = (L >= 24) ? ((size_t)1 << 27) : ((size_t)1 << 22);
  AHash H; H.init(hashSize);
  const int MAXF = 6000;
  vector<vector<uint64_t>> bucket(MAXF);
  // a state key uses the low 40 bits (mask<<8 | keeper); pack g into the high bits
  // of the bucket entry so we needn't recompute the heuristic at pop time.
  auto keyOf = [&](uint32_t mask, int keeper) -> uint64_t { return ((uint64_t)mask << 8) | (uint32_t)keeper; };
  const uint64_t KEYMASK = (1ull << 40) - 1;

  int h0 = heuristic(boxMask0);
  H.improve(keyOf(boxMask0, keeper0), 0);           // g 0
  bucket[h0].push_back(keyOf(boxMask0, keeper0));    // g packed as high bits = 0

  atomic<long long> aExpanded(0), aGenerated(0);
  atomic<int> best(INT32_MAX);
  vector<Scratch> scratch(NT);
  for (auto& sc : scratch) { sc.kdist.assign(N, 0); sc.kstamp.assign(N, 0); sc.queue.assign(N, 0); sc.kcur = 0; }
  vector<vector<pair<int, uint64_t>>> tout(NT);     // per-thread successors: (nf, entry)
  long long lastReport = 0;
  bool timedOut = false; int lbAtTimeout = 0;
  auto t0 = chrono::steady_clock::now();

  // expand one wave of same-f entries; returns the successors that land back in f.
  auto processWave = [&](vector<uint64_t>& wave, int f) -> vector<uint64_t> {
    for (auto& b : tout) b.clear();
    auto worker = [&](unsigned idx, size_t stride) {
      Scratch& sc = scratch[idx];
      auto& out = tout[idx];
      for (size_t j = idx; j < wave.size(); j += stride) {
        uint64_t entry = wave[j];
        uint64_t key = entry & KEYMASK; int g = (int)(entry >> 40);
        long vi = H.find(key);
        if (vi < 0) continue;
        uint32_t v = H.val(vi); if (v & AHash::CLOSED) continue;   // already expanded
        if ((int)(v & 0x7fffff) != g) continue;                    // stale (better g known)
        if (!H.tryClose(vi, key, (uint32_t)g)) continue;           // another thread won it
        aExpanded.fetch_add(1, memory_order_relaxed);
        uint32_t mask = (uint32_t)(key >> 8); int keeper = (int)(key & 0xff);
        keeperBFS(sc, mask, keeper);
        // all boxes placed? -> finish by walking keeper to nearest free goal
        if ((mask & goalBoxMask) == mask) {
          int fw = -1;
          for (int gc : goalCells) if (!boxAt(mask, gc) && kReach(sc, gc) && (fw < 0 || sc.kdist[gc] < fw)) fw = sc.kdist[gc];
          if (fw >= 0) { int cand = g + fw; int cb = best.load(memory_order_relaxed); while (cand < cb && !best.compare_exchange_weak(cb, cand, memory_order_relaxed)); }
        }
        // expand pushes
        uint32_t m = mask;
        while (m) {
          int b = __builtin_ctz(m); m &= m - 1; int bcell = liveCell[b];
          for (int d = 0; d < 4; d++) {
            int pf = NB[bcell * 4 + OPP[d]]; if (!kReach(sc, pf)) continue;
            int t = NB[bcell * 4 + d];
            if (t < 0 || WALL[t] || deadCell(t) || boxAt(mask, t)) continue;
            uint32_t nmask = (mask & ~(1u << b)) | (1u << liveIdx[t]);
            if (!ISGOAL[t] && isFrozen(t, nmask)) continue;
            int ng = g + sc.kdist[pf] + 1;
            uint64_t nkey = keyOf(nmask, bcell);
            if (H.improve(nkey, (uint32_t)ng)) {                   // we set a new/better g
              int nf = ng + heuristic(nmask);
              if (nf < MAXF) { out.push_back({nf, nkey | ((uint64_t)ng << 40)}); aGenerated.fetch_add(1, memory_order_relaxed); }
            }
          }
        }
      }
    };
    if (NT == 1 || wave.size() < 2000) { worker(0, 1); }           // small wave: skip thread overhead
    else { vector<thread> ths; for (unsigned t = 0; t < NT; t++) ths.emplace_back(worker, t, (size_t)NT); for (auto& th : ths) th.join(); }
    // merge per-thread successors (single-threaded, no races on bucket)
    vector<uint64_t> next;
    for (auto& buf : tout) for (auto& pr : buf) { if (pr.first == f) next.push_back(pr.second); else bucket[pr.first].push_back(pr.second); }
    return next;
  };

  for (int f = 0; f < MAXF && !timedOut; f++) {
    if (f >= best.load(memory_order_relaxed)) break;   // admissible: nothing below `best` remains
    if (bucket[f].empty()) continue;
    vector<uint64_t> wave = move(bucket[f]); bucket[f].clear();
    while (!wave.empty()) {
      wave = processWave(wave, f);
      long long ex = aExpanded.load(memory_order_relaxed);
      if (ex - lastReport >= 4000000) {
        lastReport = ex;
        double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        int bb = best.load(memory_order_relaxed);
        printf("  ...expanded %lld, f=%d gen=%lld best=%d %.0fs\n", ex, f, aGenerated.load(memory_order_relaxed), bb == INT32_MAX ? -1 : bb, el);
        fflush(stdout);
        if (el > timeout) { timedOut = true; lbAtTimeout = f; break; }
      }
    }
    bucket[f].clear(); bucket[f].shrink_to_fit();
  }

  double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
  int bestv = best.load(memory_order_relaxed);
  long long ex = aExpanded.load(memory_order_relaxed), gen = aGenerated.load(memory_order_relaxed);
  if (timedOut) { printf("RESULT %s optimal=TIMEOUT lower_bound=%d best_found=%s expanded=%lld seconds=%.1f threads=%u\n", board.c_str(), lbAtTimeout, bestv == INT32_MAX ? "none" : to_string(bestv).c_str(), ex, el, NT); return 0; }
  if (bestv != INT32_MAX) printf("RESULT %s optimal=%d expanded=%lld generated=%lld seconds=%.1f threads=%u\n", board.c_str(), bestv, ex, gen, el, NT);
  else printf("RESULT %s optimal=NONE expanded=%lld seconds=%.1f threads=%u\n", board.c_str(), ex, el, NT);
  return 0;
}
