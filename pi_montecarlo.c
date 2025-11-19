#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BASE_SEED 123456789u

typedef long long ll;

typedef struct {
    double *xs;
    double *ys;
    ll start;
    ll end;
    ll inCircle;
} WorkerLocal;

typedef struct {
    double *xs;
    double *ys;
    ll start;
    ll end;
} WorkerShared;

pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;
ll g_inCircle = 0;

double wall_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

double rand_in_minus1_1(unsigned int *seed) {
    double u = (double)rand_r(seed) / (double)RAND_MAX;
    return u * 2.0 - 1.0;
}

void generate_points(double *xs, double *ys, ll nPoints, unsigned int seed) {
    unsigned int s = seed;
    for (ll i = 0; i < nPoints; i++) {
        double x = rand_in_minus1_1(&s);
        double y = rand_in_minus1_1(&s);
        xs[i] = x;
        ys[i] = y;
    }
}

double estimate_pi_single(double *xs, double *ys, ll nPoints) {
    ll inCircle = 0;
    for (ll i = 0; i < nPoints; i++) {
        double x = xs[i];
        double y = ys[i];
        if (x * x + y * y <= 1.0) {
            inCircle++;
        }
    }
    return 4.0 * (double)inCircle / (double)nPoints;
}

void *worker_local(void *arg) {
    WorkerLocal *w = (WorkerLocal *)arg;
    ll count = 0;
    for (ll i = w->start; i < w->end; i++) {
        double x = w->xs[i];
        double y = w->ys[i];
        if (x * x + y * y <= 1.0) {
            count++;
        }
    }
    w->inCircle = count;
    return NULL;
}

double estimate_pi_multi_local(double *xs, double *ys, ll nPoints, int nThreads) {
    pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * nThreads);
    WorkerLocal *data = (WorkerLocal *)malloc(sizeof(WorkerLocal) * nThreads);

    ll base = nPoints / nThreads;
    ll extra = nPoints % nThreads;
    ll offset = 0;

    for (int i = 0; i < nThreads; i++) {
        ll len = base + (i < extra ? 1 : 0);
        data[i].xs = xs;
        data[i].ys = ys;
        data[i].start = offset;
        data[i].end = offset + len;
        data[i].inCircle = 0;
        pthread_create(&threads[i], NULL, worker_local, &data[i]);
        offset += len;
    }

    ll totalInCircle = 0;
    for (int i = 0; i < nThreads; i++) {
        pthread_join(threads[i], NULL);
        totalInCircle += data[i].inCircle;
    }

    free(threads);
    free(data);

    return 4.0 * (double)totalInCircle / (double)nPoints;
}

void *worker_shared(void *arg) {
    WorkerShared *w = (WorkerShared *)arg;
    for (ll i = w->start; i < w->end; i++) {
        double x = w->xs[i];
        double y = w->ys[i];
        if (x * x + y * y <= 1.0) {
            pthread_mutex_lock(&g_lock);
            g_inCircle++;
            pthread_mutex_unlock(&g_lock);
        }
    }
    return NULL;
}

double estimate_pi_multi_shared(double *xs, double *ys, ll nPoints, int nThreads) {
    pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * nThreads);
    WorkerShared *data = (WorkerShared *)malloc(sizeof(WorkerShared) * nThreads);

    g_inCircle = 0;

    ll base = nPoints / nThreads;
    ll extra = nPoints % nThreads;
    ll offset = 0;

    for (int i = 0; i < nThreads; i++) {
        ll len = base + (i < extra ? 1 : 0);
        data[i].xs = xs;
        data[i].ys = ys;
        data[i].start = offset;
        data[i].end = offset + len;
        pthread_create(&threads[i], NULL, worker_shared, &data[i]);
        offset += len;
    }

    for (int i = 0; i < nThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    ll totalInCircle = g_inCircle;

    free(threads);
    free(data);

    return 4.0 * (double)totalInCircle / (double)nPoints;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s <nPoints>\n", argv[0]);
        printf("  %s <nPoints> <nThreads>\n", argv[0]);
        return 1;
    }

    ll nPoints = atoll(argv[1]);

    if (nPoints <= 0) {
        printf("nPoints must be > 0\n");
        return 1;
    }

    if (argc >= 3) {
        int nThreads = atoi(argv[2]);
        if (nThreads <= 0) {
            printf("nThreads must be > 0\n");
            return 1;
        }

        double *xs = (double *)malloc(sizeof(double) * nPoints);
        double *ys = (double *)malloc(sizeof(double) * nPoints);
        if (!xs || !ys) {
            fprintf(stderr, "Cannot allocate memory\n");
            free(xs);
            free(ys);
            return 1;
        }

        generate_points(xs, ys, nPoints, BASE_SEED);

        printf("Single run mode: nPoints = %lld, nThreads = %d\n", nPoints, nThreads);

        double t1_start = wall_time();
        double pi1 = estimate_pi_single(xs, ys, nPoints);
        double t1_end = wall_time();
        double T1 = t1_end - t1_start;
        printf("\nApproach 1 (single thread)\n");
        printf("Pi ≈ %.10f, error = %.10f, time = %.6f s\n",
               pi1, fabs(pi1 - M_PI), T1);

        double t2_start = wall_time();
        double pi2 = estimate_pi_multi_local(xs, ys, nPoints, nThreads);
        double t2_end = wall_time();
        double T2 = t2_end - t2_start;
        printf("\nApproach 2 (multi-thread, local count)\n");
        printf("Pi ≈ %.10f, error = %.10f, time = %.6f s, speedup = %.3f\n",
               pi2, fabs(pi2 - M_PI), T2, T1 / T2);

        double t3_start = wall_time();
        double pi3 = estimate_pi_multi_shared(xs, ys, nPoints, nThreads);
        double t3_end = wall_time();
        double T3 = t3_end - t3_start;
        printf("\nApproach 3 (multi-thread, shared count + mutex)\n");
        printf("Pi ≈ %.10f, error = %.10f, time = %.6f s, speedup = %.3f\n",
               pi3, fabs(pi3 - M_PI), T3, T1 / T3);

        free(xs);
        free(ys);
        return 0;
    }

    printf("Benchmark mode: nPoints = %lld\n", nPoints);

    double *xs = (double *)malloc(sizeof(double) * nPoints);
    double *ys = (double *)malloc(sizeof(double) * nPoints);
    if (!xs || !ys) {
        fprintf(stderr, "Cannot allocate memory\n");
        free(xs);
        free(ys);
        return 1;
    }

    generate_points(xs, ys, nPoints, BASE_SEED);

    int REPEAT = 10;

    double baseSum = 0.0;
    double lastPiSingle = 0.0;
    for (int r = 0; r < REPEAT; r++) {
        double tStart = wall_time();
        lastPiSingle = estimate_pi_single(xs, ys, nPoints);
        double tEnd = wall_time();
        baseSum += (tEnd - tStart);
    }
    double T_single_avg = baseSum / REPEAT;
    double singleErr = fabs(lastPiSingle - M_PI);

    printf("\nApproach 1 (single thread)\n");
    printf("Pi ≈ %.10f, error = %.10f, avg time = %.6f s\n\n",
           lastPiSingle, singleErr, T_single_avg);

    int threadList[] = {2, 4, 8, 16, 32, 64, 128};
    int nConfigs = sizeof(threadList) / sizeof(threadList[0]);

    printf("=== Approach 2: multi-thread, local count ===\n");
    printf("Threads,AvgTime(s),Speedup,Pi(avg),AbsError\n");
    printf("1,%.6f,1.000,%.10f,%.10f\n",
           T_single_avg, lastPiSingle, singleErr);

    for (int i = 0; i < nConfigs; i++) {
        int N = threadList[i];
        double sum = 0.0;
        double piSum = 0.0;
        for (int r = 0; r < REPEAT; r++) {
            double tStart = wall_time();
            double pi = estimate_pi_multi_local(xs, ys, nPoints, N);
            double tEnd = wall_time();
            sum += (tEnd - tStart);
            piSum += pi;
        }
        double T_multi_avg = sum / REPEAT;
        double avgPi = piSum / REPEAT;
        double speedup = T_single_avg / T_multi_avg;
        double err = fabs(avgPi - M_PI);
        printf("%d,%.6f,%.3f,%.10f,%.10f\n",
               N, T_multi_avg, speedup, avgPi, err);
    }

    printf("\n=== Approach 3: multi-thread, shared count + mutex ===\n");
    printf("Threads,AvgTime(s),Speedup,Pi(avg),AbsError\n");
    for (int i = 0; i < nConfigs; i++) {
        int N = threadList[i];
        double sum = 0.0;
        double piSum = 0.0;
        for (int r = 0; r < REPEAT; r++) {
            double tStart = wall_time();
            double pi = estimate_pi_multi_shared(xs, ys, nPoints, N);
            double tEnd = wall_time();
            sum += (tEnd - tStart);
            piSum += pi;
        }
        double T_multi_avg = sum / REPEAT;
        double avgPi = piSum / REPEAT;
        double speedup = T_single_avg / T_multi_avg;
        double err = fabs(avgPi - M_PI);
        printf("%d,%.6f,%.3f,%.10f,%.10f\n",
               N, T_multi_avg, speedup, avgPi, err);
    }

    free(xs);
    free(ys);
    return 0;
}
