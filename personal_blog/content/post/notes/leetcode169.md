---
title: "LeetCode 169. Majority Element"
author: "Mozhang Guo"
date: 2024-02-18
lastmode: 2024-02-18
draft: false

tags: [
    "C++",
    "Boyer-Moore Voting Algorithm"
]
categories: [
    "Notes"
]

katex: true
markup: mmark

---

This problem asked us to return the majority elements that appears more than **[n/2]** times.  The problem itself is not that hard to complete with hash table or sorting algorithms.  The intuition behind this problem can be understood as to find the frequency of each element in the given array and return the most frequent element. The table below summarize the space and time complexity of these algorirthms. 

| Algorithm|Time Complexity | Space Complexity|
| :---|:----:|---:|
    | Hash Table      | O(N)       | O(N)   |
    | Sorting   | O(N)        | O(N)      |

However, the follow-up question asked us to solve the problem in linear time O(n) and constant space O(1). The Boyer-moore majority vote algorithm can be applied here to find the majority element under the constraint.

The algorithm works by first setting the majority elements and element count as 0. As we looping through the array, once we encounter the element again, we increment the counter. Otherwise, we decrement the counter. Once the counter is 0, we would set the majoirty element as the element we see now. At the end of the process, the element will be in the majorit.


```cpp
class Solution {

public:
int majorityElement(vector\<int>& nums) {
    int majnum_cnt = 0;
    int majnum_holder = 0;
    for (auto num : nums) {
    	if (majnum_cnt == 0) majnum_holder = num;
    	majnum_cnt += (num == majnum_holder) ? 1 : -1;
    }
    return majnum_holder;
}

};
```