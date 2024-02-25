---
title: "LeetCode 1188. Design Bounded Blocking Queue"
author: "Mozhang Guo"
date: 2024-02-18
lastmode: 2024-02-18
draft: false

tags: [
    "C++",
    "Concurrency"
]
categories: [
    "Notes"
]

katex: true
markup: mmark

---
# Problem Description
The problem asked us to imlement a thread safer bounded blocking queues with following methods.

* Constructor
* enqueue
* dequeue
* size

The implementation will be tested under multi-threaded condition. Each thread will either be a producer (calling enqueue) or a consumer(calling dequeue).

# Intuition
C++ offers us multiple ways of implementing the thread safe mechanism. The Conditional variable and the mutex locks can be a 
# Approach

# Complexity


# Code
Our solution using condition variable and mutex can be seen below.
```c++
class BoundedBlockingQueue {
private:
    int _capacity;
    queue<int> _queue;
    mutex _mu;
    condition_variable _cv;

public:
    BoundedBlockingQueue(int capacity) {
        _capacity = capacity;
    }
    
    void enqueue(int element) {
        std::unique_lock<mutex> lk(_mu);
        _cv.wait(lk, [this]{ return _queue.size() < _capacity;} );
        _queue.push(element);
        lk.unlock();
        _cv.notify_all();
    }
    
    int dequeue() {
        int el;
        std::unique_lock<mutex> lk(_mu);
        _cv.wait(lk, [this]{ return _queue.size() != 0; });
        el = _queue.front();
        _queue.pop();
        lk.unlock();
        _cv.notify_all();
        return el;
    }
    
    int size() {
        unique_lock<mutex> lk(_mu);
        return _queue.size();
    }
};
```