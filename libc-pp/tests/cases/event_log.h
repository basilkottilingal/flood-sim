/* The log the multi-file event tests write to. */
#ifndef EVENT_LOG_H
#define EVENT_LOG_H
#include <stdio.h>
#include <cwrap/event.h>

extern char log_ [];
extern int n_;
#define LOG(...) (n_ += snprintf (log_ + n_, 4096 - n_, __VA_ARGS__))
#endif
