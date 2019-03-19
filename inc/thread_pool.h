#ifndef _THREAD_POOL_
#define _THREAD_POOL_
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <sched.h>
#include <pthread.h>
#include <stdbool.h>

typedef void (*Fun)(void *arg);
// define a struct for passing argument to a thread function
struct Arg
{
	int coreId;
	int re_idx;
	int data_idx;
};

struct Task
{
	Fun myfun;
	void *arg;
	struct Task *next;
};

struct Thread_Pool
{
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	pthread_mutex_t mutex_flag;
	pthread_cond_t cond_flag;
	struct Task *taskHead;
	bool isClose;
	int threadNum;
	int threadNum_Idle;
	pthread_t *threadId;
};

struct pool_arg_t
{
	int coreId;
	int pool_index;
};

void pool_init(int coreId_start, int _threadNum, int pool_index);
void pool_add_task(Fun myfun, void *arg, int pool_index);
void pool_destroy(int pool_index);
void *thread_run(void *arg);

extern struct Thread_Pool *pool[12];

#endif
