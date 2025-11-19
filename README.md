Monte Carlo π Estimation with POSIX Threads

This project implements three different methods to estimate the value of π using the Monte Carlo technique. The goal is to compare the performance of a single-thread implementation with two multi-threaded variants using POSIX threads (pthread).

The program supports:

Approach 1: Single-thread baseline

Approach 2: Multi-thread with local counters (no shared variable)

Approach 3: Multi-thread with a shared global counter protected by a mutex

It can run in:

Single-run mode: run all three approaches once and print their times

Benchmark mode: sweep different thread counts and print speedup data for plotting charts

1. Monte Carlo Estimation of π

We estimate π using a classic geometric Monte Carlo method:

Consider the square 
[
−
1
,
1
]
×
[
−
1
,
1
]
[−1,1]×[−1,1] and the unit circle 
𝑥
2
+
𝑦
2
≤
1
x
2
+y
2
≤1 inscribed inside the square.

Generate nPoints random points 
(
𝑥
,
𝑦
)
(x,y) uniformly in the square.

Let inCircle be the number of points that fall inside the circle.

Because the probability of a random point lying inside the circle is:

𝑝
=
area of circle
area of square
=
𝜋
⋅
1
2
(
2
⋅
1
)
2
=
𝜋
4
,
p=
area of square
area of circle
	​

=
(2⋅1)
2
π⋅1
2
	​

=
4
π
	​

,

we can estimate π as:

𝜋
^
≈
4
×
inCircle
nPoints
.
π
^
≈4×
nPoints
inCircle
	​

.

The larger nPoints is, the closer the approximation tends to be (in expectation).

2. Program Structure

All code is contained in a single C file, for example:

pi_montecarlo.c

2.1 Common utilities

The program defines:

A fixed seed BASE_SEED = 123456789u for reproducible random points.

typedef long long ll; for 64-bit counters.

Two worker structs:

WorkerLocal for the multi-thread local-counter approach.

WorkerShared for the shared-counter approach.

Utility functions:

double wall_time() – returns wall-clock time in seconds using gettimeofday.

double rand_in_minus1_1(unsigned int *seed) – generate a random double in 
[
−
1
,
1
]
[−1,1].

void generate_points(double *xs, double *ys, ll nPoints, unsigned int seed) – fill arrays xs and ys with random points.

The same array of random points is reused by all three approaches, so they operate on identical input data.

2.2 Approach 1 – Single-thread

Function:

double estimate_pi_single(double *xs, double *ys, ll nPoints);


Behavior:

Loops through all nPoints points in arrays xs and ys.

Counts how many satisfy x*x + y*y <= 1.0.

Returns 4.0 * inCircle / nPoints.

In benchmark mode, this function is executed several times to compute an average runtime T_single_avg. This average is used as the baseline for computing speedup in the multi-thread approaches.

2.3 Approach 2 – Multithread with Local Counters

Functions:

void *worker_local(void *arg);
double estimate_pi_multi_local(double *xs, double *ys, ll nPoints, int nThreads);


Behavior:

estimate_pi_multi_local splits the index range [0, nPoints) into nThreads contiguous segments.

Each thread gets a WorkerLocal struct with:

xs, ys pointers.

start and end indices.

A local counter inCircle.

worker_local:

Loops from start to end,

Checks x*x + y*y <= 1.0,

Increments the local inCircle only.

After joining all threads, the main thread sums up all inCircle values and computes π.

Key point:

There is no shared counter and no lock in the inner loop.

Threads are completely independent during counting.

2.4 Approach 3 – Multithread with Shared Global Counter + Mutex

Functions:

void *worker_shared(void *arg);
double estimate_pi_multi_shared(double *xs, double *ys, ll nPoints, int nThreads);


Behavior:

Similar to Approach 2 in how the points are divided into segments.

However, each thread increments a shared global counter g_inCircle.

To avoid race conditions, each increment is protected by a mutex g_lock:

pthread_mutex_lock(&g_lock);
g_inCircle++;
pthread_mutex_unlock(&g_lock);


After joining all threads, g_inCircle holds the total number of points inside the circle, and π is computed from this global counter.

This approach is correct (no data races), but introduces heavy synchronization overhead, because every hit inside the circle involves a lock/unlock pair.

3. Building the Program

Requirements:

A C compiler (e.g. gcc or clang)

POSIX threads (pthread)

Math library (libm)

Typical build command:

gcc pi_montecarlo.c -o pi_montecarlo -lpthread -lm -O2


Options:

-lpthread – link against the pthread library.

-lm – link against the math library (for M_PI, fabs, etc.).

-O2 – enable compiler optimizations.

4. Running the Program

The executable supports two modes, depending on the number of arguments.

4.1 Single-run mode (3 approaches once)

Usage:

./pi_montecarlo <nPoints> <nThreads>


Example:

./pi_montecarlo 10000000 8


This will:

Generate nPoints random points once.

Run:

Approach 1 (single thread),

Approach 2 (multi-thread, local counters) with nThreads,

Approach 3 (multi-thread, shared counter + mutex) with nThreads.

Print:

Estimated π,

Absolute error |π_est − π_true|,

Execution time for each approach,

Speedup of each multi-thread approach relative to the single-thread baseline.

4.2 Benchmark mode (data for speedup charts)

Usage:

./pi_montecarlo <nPoints>


Example:

./pi_montecarlo 10000000


This will:

Generate nPoints random points once.

Benchmark Approach 1 (single-thread) multiple times to get T_single_avg.

Benchmark Approach 2 (local counters) for a list of thread counts, e.g. {1, 2, 4, 8, 16, 32, 64, 128}.

Benchmark Approach 3 (shared global counter) for {2, 4, 8, 16, 32, 64, 128}.

Print tables in CSV-like format:

For Approach 1:

Approach 1 (single thread)
Pi ≈ 3.1418992000, error = 0.0003065464, avg time = 0.000930 s


For Approach 2:

=== Approach 2: multi-thread, local count ===
Threads,AvgTime(s),Speedup,Pi(avg),AbsError
1,0.000930,1.000,3.1418992000,0.0003065464
2,0.002729,0.341,3.1418992000,0.0003065464
4,0.002601,0.358,3.1418992000,0.0003065464
...


For Approach 3:

=== Approach 3: multi-thread, shared count + mutex ===
Threads,AvgTime(s),Speedup,Pi(avg),AbsError
2,0.099510,0.009,3.1418992000,0.0003065464
4,0.231788,0.004,3.1418992000,0.0003065464
...


You can copy these tables directly into Excel / Google Sheets to draw speedup charts.

5. Plotting Speedup Charts

To plot a speedup curve similar to the assignment:

Copy the output table for Approach 2 (or Approach 3) into a CSV or spreadsheet.

Use columns:

X axis: Threads (e.g. 1, 2, 4, 8, 16, 32, 64, 128)

Y axis: Speedup

Use a line chart (optionally with markers).
In some plots, the X axis can be scaled logarithmically (base 2) to match the visual style in “Exploring Multithreaded Performance”.

You will typically observe:

For Approach 2:

Speedup values < 1 (e.g. 0.30–0.36), meaning the multi-thread version is slower than the single-thread baseline due to thread-creation and scheduling overhead.

For Approach 3:

Speedup values around 0.004–0.009, meaning the shared-counter version with a mutex is much slower because every increment requires acquiring and releasing a lock.

6. Interpretation and Discussion (Short Summary)

All three approaches produce the same π approximation, because they use the same set of random points and perform the same geometric test.

Approach 1 is simple and efficient for moderate nPoints.

Approach 2 eliminates data races by using per-thread local counters, but does not gain speedup on this workload: the computational work per point is too small compared to the thread-management overhead.

Approach 3 demonstrates a correct but inefficient shared-variable design: updating a global counter with a mutex in every iteration causes severe contention and effectively serialises the computation.

A typical optimisation for Approach 3 is to use local accumulation plus one mutex update per thread: each thread keeps its own local count and only acquires the lock once at the end to add to the global counter. This preserves the shared variable g_inCircle while drastically reducing synchronization overhead.
