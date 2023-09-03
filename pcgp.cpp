#ifdef _WIN32
#define FSEEK _fseeki64
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#else
#define _FILE_OFFSET_BITS 64
#define FSEEK fseek
#endif

#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <climits>
#include <cfloat>
#include <cctype>

#include <omp.h>

#include "pcgpmem.h"

#define LINE_MAX 4096
#define MAX_THREADS 12

void die(const char *msg) {
    if (errno) {
        perror(msg);
    } else if (msg) {
        fprintf(stderr, "%s.\n", msg);
    }
    exit(EXIT_FAILURE);
}

bool str2int(int *out, const char *str, int base) {
    char *end = nullptr;
    if (str[0] == '\0' || isspace(str[0]))
        return false;
    errno = 0;
    long l = strtol(str, &end, base);
    if (l > INT_MAX || l < INT_MIN
        || (errno == ERANGE && (l == LONG_MAX || l == LONG_MIN))
        || *end != '\0') {
        return false;
    } else {
        *out = static_cast<int>(l);
        return true;
    }
}

typedef struct {
    int n;      // order
    int k;      // number of jumps
    int so;     // number of static jumps
    int *s;     // jumps
} Graph;

typedef struct {
    int     diam;
    double  aspl;
} GraphProp;

typedef struct {
    int u;
    int v;
} VertPair;

void adjVert(VertPair *__restrict p, const Graph *__restrict g, int u, int i) {
    p->u = (u + g->s[i]) % g->n;
    p->v = (u - g->s[i] + g->n) % g->n;
}

void circulantBFS(Arena *__restrict arena, GraphProp *__restrict prop, const Graph *__restrict g) {
    int *__restrict dist = static_cast<int*>(arenaAlloc(arena, g->n * sizeof(*dist)));
    int *__restrict queue = static_cast<int*>(arenaAlloc(arena, g->n * sizeof(*queue)));
    dist[0] = 0;
    for (int i = 1; i < g->n; ++i) {
        dist[i] = INT_MAX;
    }
    queue[0] = 0;
    int qr = 0;
    int qw = 1;
    int vc = 1;
    while (vc < g->n && qr != qw) {
        int u = queue[++qr];
        int d = dist[u] + 1;
        for (int i = 0; i < g->k; ++i) {
            VertPair uv;
            adjVert(&uv, g, u, i);
            if (dist[uv.u] == INT_MAX) {
                queue[++qw] = uv.u;
                dist[uv.u] = dist[(g->n - uv.u) % g->n] = d;
                vc += 2;
            }
            if (dist[uv.v] == INT_MAX) {
                queue[++qw] = uv.v;
                dist[uv.v] = dist[(g->n - uv.v) % g->n] = d;
                vc += 2;
            }
        }
    }
    int diam = 0;
    int dist_sum = 0;
    for (int i = 0; i < g->n; ++i) {
        diam = diam > dist[i] ? diam : dist[i];
        dist_sum += dist[i];
    }
    prop->diam = diam;
    prop->aspl = (double)dist_sum / (g->n - 1);
}

typedef struct {
    int a, b;
    int g;
} KLPair;

int KernighanLinC(const Graph *graph, const int *part, unsigned a, unsigned b) {
    for (int i = 0; i < graph->k; ++i) {
        VertPair uv;
        adjVert(&uv, graph, static_cast<int>(a), i);
        if (uv.v == static_cast<int>(b) || uv.u == static_cast<int>(b)) {
            return part[a] != part[b] ? 1 : -1;
        }
    }
    return 0;
}

int KernighanLinPartitionCost(const Graph *graph, const int *part) {
    int cost = 0;
    for (int u = 0; u < graph->n; ++u) {
        for (int i = 0; i < graph->k; ++i) {
            const int v = (u + graph->s[i]) % graph->n;
            if (part[u] != part[v]) {
                ++cost;
            }
        }
    }
    return cost;
}

int circulantKernighanLin(Arena *__restrict arena, const Graph *__restrict graph) {
    int *__restrict V = static_cast<int*>(arenaAlloc(arena, graph->n * sizeof(*V)));
    int *__restrict P = static_cast<int*>(arenaAlloc(arena, graph->n * sizeof(*P)));
    int *__restrict D = static_cast<int*>(arenaAlloc(arena, graph->n * sizeof(*D)));
    int *__restrict G_sum = static_cast<int*>(arenaAlloc(arena, graph->n * sizeof(*G_sum)));
    int G_sum_max = 0;
    int GAB_size = 0;
    KLPair *__restrict GAB = static_cast<KLPair*>(arenaAlloc(arena, graph->n * sizeof(*GAB)));
    for (int i = 0; i < graph->n / 2; ++i) {
        P[i] = 0;
    }
    for (int i = graph->n / 2; i < graph->n; ++i) {
        P[i] = 1;
    }
    do {
        GAB_size = 0;
        for (int u = 0; u < graph->n; ++u) { //todo
            V[u] = 0;
            D[u] = 0;

            #pragma omp parallel for default(none) shared(graph, u, P) reduction(+:D[u])
            for (int i = 0; i < graph->k; ++i) {
                VertPair uv;
                adjVert(&uv, graph, u, i);
                D[u] += (P[u] != P[uv.u] ? 1 : -1);
                if (uv.u != uv.v) {
                    D[u] += (P[u] != P[uv.v] ? 1 : -1);
                }
            }
        }
        {
            for (int i = 0; i < graph->n / 2; ++i) {
                int g_max[3];
                g_max[0] = INT_MIN;
                g_max[1] = g_max[2] = 0;

                #pragma omp parallel for default(none) shared(graph, V, P, D, g_max)
                for (int a = 0; a < graph->n; ++a) { //todo
                    if (!V[a] && (P[a] == 0)) {
                        for (int b = 0; b < graph->n; ++b) {
                            if (!V[b] && (P[b] == 1)) {
                                int g = D[a] + D[b] - 2 * KernighanLinC(graph, P, a, b);
                                if (g > g_max[0]) {
                                    g_max[0] = g;
                                    g_max[1] = a;
                                    g_max[2] = b;
                                }
                            }
                        }
                    }
                }

                #pragma omp master
                V[g_max[1]] = V[g_max[2]] = 1;
                KLPair *gab = &GAB[++GAB_size];
                gab->g = g_max[0];
                gab->a = g_max[1];
                gab->b = g_max[2];

                #pragma omp parallel for default(none) shared(graph, V, P, D, g_max)
                for (int u = 0; u < graph->n; ++u) { //todo
                    if (V[u]) continue;
                    else if (P[u]) {
                        int c_yb = KernighanLinC(graph, P, u, g_max[2]);
                        int c_ab = KernighanLinC(graph, P, u, g_max[1]);
                        D[u] += 2 * c_yb - 2 * c_ab;
                    } else {
                        int c_xa = KernighanLinC(graph, P, u, g_max[1]);
                        int c_xb = KernighanLinC(graph, P, u, g_max[2]);
                        D[u] += 2 * c_xa - 2 * c_xb;
                    }
                }
            }
        }
        int k = 0;
        G_sum_max = G_sum[0] = GAB[0].g;
        for (int i = 1; i < GAB_size; ++i) {
            G_sum[i] = G_sum[i - 1] + GAB[i].g;
            if (G_sum[i] > G_sum_max) {
                G_sum_max = G_sum[i];
                k = i;
            }
        }
        if (G_sum_max > 0) {
            for (int i = 0; i <= k; ++i) {
                KLPair *uv = &GAB[i];
                int tmp = P[uv->a];
                P[uv->a] = P[uv->b];
                P[uv->b] = tmp;
            }
        }
    } while (G_sum_max > 0);
    return KernighanLinPartitionCost(graph, P);
}

enum PcgpStage {
    PCGP_STAGE_BFS,
    PCGP_STAGE_KL,
};

typedef struct {
    int                 stage;
    int                 best_diam;
    double              best_aspl;
    int                 best_bisect_cost;
    unsigned long long  graph_count;
    unsigned long long  graph_count_max;
    unsigned long long  graph_count_best;
    size_t              bfs_log_count;
    size_t              kl_log_count;
} ScanState;

int stateFileWrite(ScanState *scan, Graph *g, FILE *file) {
    void *buf = g->s + g->so;
    size_t len = g->k - g->so;
    if (FSEEK(file, 0L, SEEK_SET)) {
        return 1;
    } else if (fwrite(&scan, sizeof(scan), 1, file) < 1) {
        return 1;
    } else if (fwrite(buf, sizeof(*g->s), len, file) < len) {
        return 1;
    } else {
        return 0;
    }
}

#define STATE_UPDATE_TIME 5.0

enum PcgpMode {
    PCGP_MODE_SCAN,
    PCGP_MODE_IMMEDIATE,
};

bool graphCheck(const Graph *g) {
    if (g->n < 0 || g->k <= 0 || g->k > g->n / 2 || g->so < 0 || g->so > g->k) {
        return true;
    }
    return false;
}

int main(int argc, char **argv) {
    int mode = -1;
    if (argc >= 1) {
        if (strcmp(argv[1], "s") == 0) {
            mode = PCGP_MODE_SCAN;
        } else if (strcmp(argv[1], "i") == 0) {
            mode = PCGP_MODE_IMMEDIATE;
        }
    }
    --argc;
    ++argv;
    omp_set_num_threads(MAX_THREADS);
    switch (mode) {
    case PCGP_MODE_SCAN: {
        Graph g;
        if (argc != 4 || !(str2int(&g.n, argv[1], 10) && str2int(&g.k, argv[2], 10) && str2int(&g.so, argv[3], 10)))
            die("Invalid arguments");
        if (graphCheck(&g))
            die("Invalid argument values");
        g.s = static_cast<int*>(malloc(g.k * sizeof(*g.s)));
        if (!g.s)
            goto error;
        const size_t GS_SIZE = (g.k - g.so) * sizeof(*g.s);
        ScanState scan;
        scan.stage = PCGP_STAGE_BFS;
        scan.best_diam = INT_MAX;
        scan.best_aspl = DBL_MAX;
        scan.best_bisect_cost = 0;
        unsigned long long graph_count_prev = 0;
        FILE *state_file = nullptr;
        FILE *bfs_file = nullptr;
        FILE *kl_file = nullptr;
        {/* restore scan state from file or initialize */
            bool restored = false;
            char bfs_file_name[0x100];
            {
                int len = snprintf(bfs_file_name, sizeof(bfs_file_name), "%d-%d-%d.bfs.bin", g.n, g.k, g.so);
                if (len < 0 || (unsigned)len >= sizeof(bfs_file_name))
                    goto error;
            }
            char kl_file_name[0x100];
            {
                int len = snprintf(kl_file_name, sizeof(kl_file_name), "%d-%d-%d.kl.bin", g.n, g.k, g.so);
                if (len < 0 || (unsigned)len >= sizeof(kl_file_name))
                    goto error;
            }
            char state_file_name[0x100];
            {
                int len = snprintf(state_file_name, sizeof(state_file_name), "%d-%d-%d.state.bin", g.n, g.k, g.so);
                if (len < 0 || (unsigned)len >= sizeof(state_file_name))
                    goto error;
            }
            state_file = fopen(state_file_name, "rb+");
            if (state_file) {
                if (fread(&scan, sizeof(scan), 1, state_file) == 1) {
                    for (int i = 0; i < g.so; ++i) {
                        g.s[i] = i + 1;
                    }
                    void *buf = g.s + g.so;
                    size_t len = g.k - g.so;
                    if (fread(buf, sizeof(*g.s), len, state_file) == len) {
                        restored = true;
                    }
                }
            }
            if (restored) {
                bfs_file = fopen(bfs_file_name, "rb+");
                if (!bfs_file)
                    goto error;
                if (FSEEK(bfs_file, static_cast<long>(scan.bfs_log_count * GS_SIZE), SEEK_SET))
                    goto error;
                kl_file = fopen(kl_file_name, "rb+");
                if (!kl_file)
                    goto error;
                if (FSEEK(kl_file, static_cast<long>(scan.kl_log_count * GS_SIZE), SEEK_SET))
                    goto error;
                graph_count_prev = scan.graph_count;
            } else {/* initialize scan state */
                bfs_file = fopen(bfs_file_name, "wb+");
                if (!bfs_file)
                    goto error;
                kl_file = fopen(kl_file_name, "wb+");
                if (!kl_file)
                    goto error;
                state_file = fopen(state_file_name, "wb+");
                if (!state_file)
                    goto error;
                {/* compute max graph count as C(n,k) (overflow not tested) */
                    int n = g.n / 2 - g.so;
                    int k = g.k - g.so;
                    if (k * 2 > n) k = n - k;
                    if (k > n) {
                        scan.graph_count_max = 0;
                    } else if (k == 0) {
                        scan.graph_count_max = 1;
                    } else {
                        long long c = n;
                        for (long long i = 2; i <= k; ++i) {
                            c *= (n - i + 1);
                            c /= i;
                        }
                        scan.graph_count_max = c;
                    }
                }
                for (int i = 0; i < g.k; ++i) {
                    g.s[i] = i + 1;
                }
                if (stateFileWrite(&scan, &g, state_file))
                    goto error;
            }
            fprintf(stderr, "N=%d K=%d C=%d (%llu graphs)\n", g.n, g.k, g.so, scan.graph_count_max);
            if (restored) {
                fprintf(stderr, "Restored scan state at stage %d graph %llu\n", scan.stage, scan.graph_count);
            }
        }
        FILE *log_file = nullptr;
        {
            char file_name[0x100];
            int len = snprintf(file_name, sizeof(file_name), "%d-%d-%d.log.csv", g.n, g.k, g.so);
            if (len < 0 || (unsigned)len >= sizeof(file_name))
                goto error;
            log_file = fopen(file_name, "w");
            if (!log_file)
                goto error;
            if (fputs("graph_count,elapsed,eta,best_diam,best_aspl\n", log_file) == EOF)
                goto error;
        }
        Arena arena;
        {
            size_t arena_size = 64 * g.n * sizeof(*g.s);
            void *arena_mem = malloc(arena_size);
            if (!arena_mem)
                goto error;
            arenaInit(&arena, arena_mem, arena_size);
        }
        clock_t clock_start = 0;
        clock_t clock_prev = 0;
        switch (scan.stage) {
        case PCGP_STAGE_BFS: {
            clock_start = clock_prev = clock();
            bool running = true;
            while (running) {
                {/* compute and write values for S */
                    GraphProp prop;
                    circulantBFS(&arena, &prop, &g);
                    if (prop.diam <= scan.best_diam && prop.aspl <= scan.best_aspl) {
                        if (prop.aspl < scan.best_aspl) {
                            scan.best_aspl = prop.aspl;
                            scan.best_diam = prop.diam;
                            scan.graph_count_best = 1;
                        } else {
                            ++scan.graph_count_best;
                        }
                        ++scan.bfs_log_count;
                        {/* write only the variable part of S */
                            void *buf = g.s + g.so;
                            size_t len = g.k - g.so;
                            if (fwrite(buf, sizeof(*g.s), len, bfs_file) < len)
                                goto error;
                        }
                    }
                }
                {/* compute the next lexicographical S or stop */
                    int n = g.n / 2;
                    int k = g.k - g.so;
                    int *s = g.s + g.so;
                    bool changed = false;
                    for (int i = 1; i <= k; ++i) {
                        if (s[k - i] <= n - i) {
                            ++s[k - i];
                            for (int j = k - i + 1; j < k; ++j) {
                                s[j] = s[j - 1] + 1;
                            }
                            changed = true;
                            break;
                        }
                    }
                    running = changed;
                }
                ++scan.graph_count;
                clock_t clock_now = clock();
                double elapsed = static_cast<double>(clock_now - clock_prev) / CLOCKS_PER_SEC;
                if (!running || elapsed >= STATE_UPDATE_TIME) {
                    {/* print logs */
                        clock_prev = clock_now;
                        unsigned long long graph_count_rel = scan.graph_count - graph_count_prev;
                        graph_count_prev     = scan.graph_count;
                        double elapsed_total = (double)(clock_now - clock_start) / CLOCKS_PER_SEC;
                        double gps           = elapsed != 0 ? (double)graph_count_rel / elapsed : (double)graph_count_rel;
                        double eta           = (double)(scan.graph_count_max - scan.graph_count) / gps;
                        double scan_prc      = (double)scan.graph_count / (double)scan.graph_count_max * 100.0;
                        int rc = 0;
                        rc = fprintf(stderr, "[Stage 0] %.2f%% (%lld/%lld) ETA=%.fs Elapsed=%gs GPS=%g Diam=%d ASPL=%g\n",
                                     scan_prc, scan.graph_count, scan.graph_count_max,
                                     eta, elapsed_total, gps, scan.best_diam, scan.best_aspl);
                        if (rc < 0)
                            goto error;
                        rc = fprintf(log_file, "%llu,%.f,%.f,%d,%g\n",
                                     scan.graph_count, elapsed_total, eta, scan.best_diam, scan.best_aspl);
                        if (rc < 0)
                            goto error;
                    }
                    {/* checkpoint the state file */
                        if (FSEEK(state_file, 0L, SEEK_SET))
                            goto error;
                        if (fwrite(&scan, sizeof(scan), 1, state_file) < 1)
                            goto error;
                        {
                            void *buf = g.s + g.so;
                            size_t len = g.k - g.so;
                            if (fwrite(buf, sizeof(*g.s), len, state_file) < len)
                                goto error;
                        }
                    }
                }
                arenaFree(&arena);
            }
            {/* write BFS pass results */
                FILE *output_file = nullptr;
                {
                    char file_name[0x100];
                    int len = snprintf(file_name, sizeof(file_name), "%d-%d-%d.output-stage0.txt", g.n, g.k, g.so);
                    if (len < 0 || (unsigned)len >= sizeof(file_name))
                        goto error;
                    output_file = fopen(file_name, "w");
                    if (!output_file)
                        goto error;
                    fprintf(stderr, "Saving stage 0 output to \"%s\"\n", file_name);
                }
                fprintf(output_file, "# Version PCGP 1.0\n");
                fprintf(output_file, "# Stage 0\n");
                fprintf(output_file, "# Parameters N=%d K=%d C=%d\n", g.n, g.k, g.so);
                fprintf(output_file, "# Values ASPL=%g Diameter=%d\n", scan.best_aspl, scan.best_diam);
                if (FSEEK(bfs_file, static_cast<long>((scan.bfs_log_count - scan.graph_count_best) * GS_SIZE), SEEK_SET))
                    goto error;
                {
                    void *buf = g.s + g.so;
                    size_t len = g.k - g.so;
                    size_t ele_read = 0;
                    while ((ele_read = fread(buf, sizeof(*g.s), len, bfs_file))) {
                        if (ele_read == len) {
                            for (int i = 0; i < g.k - 1; ++i) {
                                fprintf(output_file, "%d ", g.s[i]);
                            }
                            fprintf(output_file, "%d\n", g.s[g.k - 1]);
                        } else {
                            goto error;
                        }
                    }
                    if (fwrite(buf, sizeof(*g.s), len, state_file) < len)
                        goto error;
                }
            }
            scan.stage = PCGP_STAGE_KL;
            scan.graph_count = 0;
            scan.graph_count_max = scan.graph_count_best;
            scan.graph_count_best = 0;
            graph_count_prev = 0;
            if (FSEEK(state_file, 0L, SEEK_SET))
                goto error;
            if (fwrite(&scan, sizeof(scan), 1, state_file) < 1)
                goto error;
        }
        case PCGP_STAGE_KL: {
            if (FSEEK(bfs_file, static_cast<long>((scan.bfs_log_count - scan.graph_count_max + scan.graph_count) * GS_SIZE), SEEK_SET))
                goto error;
            {
                void *buf = g.s + g.so;
                size_t len = g.k - g.so;
                size_t ele_read = 0;
                clock_start = clock_prev = clock();
                while ((ele_read = fread(buf, sizeof(*g.s), len, bfs_file))) {
                    if (ele_read == len) {
                        int bisect_cost = circulantKernighanLin(&arena, &g);
                        if (bisect_cost >= scan.best_bisect_cost) {
                            if (bisect_cost > scan.best_bisect_cost) {
                                scan.best_bisect_cost = bisect_cost;
                                scan.graph_count_best = 1;
                            } else {
                                ++scan.graph_count_best;
                            }
                            ++scan.kl_log_count;
                            {
                                buf = g.s + g.so;
                                len = g.k - g.so;
                                if (fwrite(buf, sizeof(*g.s), len, kl_file) < len)
                                    goto error;
                            }
                        }
                        ++scan.graph_count;
                        clock_t clock_now = clock();
                        double elapsed = static_cast<double>(clock_now - clock_prev) / CLOCKS_PER_SEC;
                        if ((scan.graph_count == scan.graph_count_max) || elapsed >= STATE_UPDATE_TIME) {
                            {/* print logs */
                                clock_prev = clock_now;
                                unsigned long long graph_count_rel = scan.graph_count - graph_count_prev;
                                graph_count_prev     = scan.graph_count;
                                double elapsed_total = (double)(clock_now - clock_start) / CLOCKS_PER_SEC;
                                double gps           = (double)graph_count_rel / elapsed;
                                double eta           = (double)(scan.graph_count_max - scan.graph_count) / gps;
                                double scan_prc      = (double)scan.graph_count / (double)scan.graph_count_max * 100.0;
                                int rc = 0;
                                rc = fprintf(stderr, "[Stage 1] %.2f%% (%lld/%lld) ETA=%.fs Elapsed=%gs GPS=%g BISECT_COST=%d\n",
                                             scan_prc, scan.graph_count, scan.graph_count_max,
                                             eta, elapsed_total, gps, scan.best_bisect_cost);
                                if (rc < 0)
                                    goto error;
                            }
                            {/* checkpoint the state file */
                                if (FSEEK(state_file, 0L, SEEK_SET))
                                    goto error;
                                if (fwrite(&scan, sizeof(scan), 1, state_file) < 1)
                                    goto error;
                            }
                        }
                    } else {
                        goto error;
                    }
                    arenaFree(&arena);
                }
            }
            {/* write KL pass results */
                FILE *output_file = nullptr;
                {
                    char file_name[0x100];
                    int len = snprintf(file_name, sizeof(file_name), "%d-%d-%d.output-stage1.txt", g.n, g.k, g.so);
                    if (len < 0 || (unsigned)len >= sizeof(file_name))
                        goto error;
                    output_file = fopen(file_name, "w");
                    if (!output_file)
                        goto error;
                    fprintf(stderr, "Saving stage 1 output to \"%s\"\n", file_name);
                }
                fprintf(output_file, "# Version PCGP 1.0\n");
                fprintf(output_file, "# Stage 1\n");
                fprintf(output_file, "# Parameters N=%d K=%d C=%d\n", g.n, g.k, g.so);
                fprintf(output_file, "# Values ASPL=%g Diameter=%d BisectCost=%d\n", scan.best_aspl, scan.best_diam, scan.best_bisect_cost);
                if (FSEEK(kl_file, static_cast<long>((scan.kl_log_count - scan.graph_count_best) * GS_SIZE), SEEK_SET))
                    goto error;
                {
                    void *buf = g.s + g.so;
                    size_t len = g.k - g.so;
                    size_t ele_read = 0;
                    while ((ele_read = fread(buf, sizeof(*g.s), len, kl_file))) {
                        if (ele_read == len) {
                            for (int i = 0; i < g.k - 1; i++) {
                                fprintf(output_file, "%d ", g.s[i]);
                            }
                            fprintf(output_file, "%d\n", g.s[g.k - 1]);
                        } else {
                            goto error;
                        }
                    }
                }
            }
        }
        }
        fprintf(stderr, "Done\n");
    } break;
    case PCGP_MODE_IMMEDIATE: {
        Arena arena;
        size_t arena_size = 1 << 26;
        {
            void *arena_mem = malloc(arena_size);
            if (!arena_mem)
                goto error;
            arenaInit(&arena, arena_mem, arena_size);
        }
        int graph_buf[BUFSIZ];
        char line_buf[LINE_MAX];
        while (fgets(line_buf, LINE_MAX, stdin)) {
            if (ferror(stdin))
                goto error;
            Graph g;
            {
                int n, len, i = 0;
                char *line = line_buf;
                while (sscanf(line, "%d%n", &n, &len) == 1) {
                    graph_buf[i] = n;
                    line += len;
                    ++i;
                }
                g.n = graph_buf[0];
                g.k = i - 1;
                g.s = graph_buf + 1;
            }
            if (graphCheck(&g)) {
                fprintf(stdout, "Invalid circulant parameters\n");
                continue;
            }
            {
                GraphProp prop;
                circulantBFS(&arena, &prop, &g);
                int bisect_cost = circulantKernighanLin(&arena, &g);
                fprintf(stdout, "%d %g %d\n", prop.diam, prop.aspl, bisect_cost);
            }
            arenaFree(&arena);
        }
    } break;
    default:
        die("Invalid mode argument");
    }
    return EXIT_SUCCESS;
error:
    die(nullptr);
}
