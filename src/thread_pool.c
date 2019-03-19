#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include "thread_pool.h"

struct Thread_Pool *pool[12];
static int coreId[72];

static int stick_this_thread_to_core(int core_id)
{
	int num_cores = sysconf(_SC_NPROCESSORS_ONLN); //read cpu_core number
	if (core_id < 0 || core_id >= num_cores)
		return EINVAL;

	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(core_id, &cpuset);

	pthread_t current_thread = pthread_self();
	return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

void pool_init(int coreId_start, int _threadNum, int pool_index)
{
	int j = pool_index;
	struct pool_arg_t *pool_arg = (struct pool_arg_t *)malloc(sizeof(struct pool_arg_t) * 72);

	pool[j] = (struct Thread_Pool *)malloc(sizeof(struct Thread_Pool));
	assert(pool[j] != NULL);

	pthread_mutex_init(&(pool[j]->mutex), NULL);
	pthread_cond_init(&(pool[j]->cond), NULL);
	pool[j]->taskHead = NULL;
	pool[j]->isClose = false;
	pool[j]->threadNum = _threadNum;
	pool[j]->threadId = (pthread_t *)malloc(sizeof(pthread_t) * pool[j]->threadNum);

	int i;
	//int coreId_start = 0;
	for (i = 0; i < pool[j]->threadNum; ++i)
	{
		coreId[i] = coreId_start + i;
		// printf("coreId_start=%d coreId=%d\n", coreId_start, coreId[i]);
		pool_arg[i].coreId = coreId[i];
		pool_arg[i].pool_index = pool_index;
		if (pthread_create(&(pool[j]->threadId[i]), NULL, thread_run, (void *)&pool_arg[i]))
		{
			printf("pthread_creat failed!\n");
			return;
		}
	}
}

void pool_add_task(Fun _myfun, void *_arg, int pool_index)
{
	int j = pool_index;
	//构造一个新任务
	struct Task *newTask = (struct Task *)malloc(sizeof(struct Task));
	newTask->myfun = _myfun;
	newTask->arg = _arg;
	newTask->next = NULL; //别忘置空
	struct Task *head;

	//将任务加到任务链表中
	pthread_mutex_lock(&(pool[j]->mutex));
	head = pool[j]->taskHead;
	if (head == NULL)
		pool[j]->taskHead = newTask;
	else
	{
		while (head->next)
			head = head->next;
		head->next = newTask;
	}
	//printf("newTask: %d\n", pthread_self());
	pthread_mutex_unlock(&(pool[j]->mutex));
	pthread_cond_signal(&(pool[j]->cond));
}

void pool_destroy(int pool_index)
{
	int j = pool_index;
	if (pool[j]->isClose == true) //防止多次调用该函数
		return;
	pool[j]->isClose = true;
	//唤醒所有等待线程，然后销毁线程池
	pthread_cond_broadcast(&(pool[j]->cond));

	//回收线程
	int i;
	for (i = 0; i < pool[j]->threadNum; ++i)
		pthread_join(pool[j]->threadId[i], NULL);
	free(pool[j]->threadId);

	//销毁任务链表
	struct Task *tmpTask;
	while (pool[j]->taskHead != NULL)
	{
		tmpTask = pool[j]->taskHead;
		pool[j]->taskHead = pool[j]->taskHead->next;
		free(tmpTask);
	}

	//销毁条件变量与互斥量
	pthread_mutex_destroy(&(pool[j]->mutex));
	pthread_cond_destroy(&(pool[j]->cond));

	free(pool[j]);
	//释放内存后将指针置空
	pool[j] = NULL;
}

void *thread_run(void *_arg)
{
	//printf("thread %d is ready\n", pthread_self());
	struct Task *curTask;
	struct pool_arg_t pool_arg = *((struct pool_arg_t *)_arg);
	int coreId = pool_arg.coreId;
	int j = pool_arg.pool_index;
	//printf("I use core %d\n", coreId);
	//stick the thread to coreId
	// printf("coreId = %d\n", coreId);
	if (stick_this_thread_to_core(coreId))
		printf("Stick to core %d is failed!\n", coreId);

	while (1)
	{
		pthread_mutex_lock(&(pool[j]->mutex));
		while (pool[j]->taskHead == NULL && pool[j]->isClose == false)
		{
			//printf("thread %d is waiting\n", pthread_self());
			pthread_cond_wait(&(pool[j]->cond), &(pool[j]->mutex));
			//printf("thread: %d wakes up, taskHead: %d\n", pthread_self(), pool[j] -> taskHead);
		}
		if (pool[j]->taskHead == NULL && pool[j]->isClose == true) //销毁线程池时保证任务链表已空
		{
			pthread_cond_broadcast(&(pool[j]->cond));
			pthread_mutex_unlock(&(pool[j]->mutex));
			//printf("thread %d is over\n", pthread_self());
			pthread_exit(NULL);
		}
		//printf("thread %d is going to work\n", pthread_self());
		//printf("thread: %d wakes up, taskHead: %d\n", pthread_self(), pool[j] -> taskHead);
		// assert(pool[j]->taskHead != NULL);

		curTask = pool[j]->taskHead;
		pool[j]->taskHead = pool[j]->taskHead->next;
		pthread_mutex_unlock(&(pool[j]->mutex));
		//执行任务函数
		(curTask->myfun)(curTask->arg);
		free(curTask);
		curTask = NULL;
	}
}
