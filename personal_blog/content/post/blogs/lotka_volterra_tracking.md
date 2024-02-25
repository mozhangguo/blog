---
title: "Lotka Volterra System Tracking with Kalman Filter"
author: "Mozhang Guo"
date: 2024-02-17
lastmode: 2024-02-18
draft: true

tags: [
    "MATLAB",
    "EKF",
    "KF",
    "UKF",
    "DSP",
    "Lotka Volterrra System"
]
categories: [
    "Blog"
]

katex: true
markup: mmark

---

# Introduction
The Kalman Filter provides a reliable solution in system tracking for linear systems. It estimates the system states and make predictions about unknown variables with noised measurements. However, for non-linear systems, the Kalman filter does not perform well since the system cannot be expressed in linear matrix form. Other Kalman Filter based sub-variant algorithms are designed for these non-linear Gaussian problems, such as Linearized Kalman Filter, Extended Kalman Filter, and Unscented Kalman Filter. In this article, our goal will be to implement and compare these algorithms for the non-linear system.

Although both Extended Kalman Filter and Unscented Kalman Filter perform well on non-linear system tracking, there are still debates on which one to choose in practical applications. Depending on the system complexity and algorithm design, the complexity of UKF and EKF are different. In the following section, we have reviewed some practical applications of UKF and EKF. These studies provide a comparative analysis of different Kalman Filter based algorithms under different physical scenarios.

# System Model

## Non-linear System

We consider a system described as classical predator-prey population dynamic system. It is considered in \cite{kovvali2013introduction},\cite{chow2004fitting},\cite{chow2007unscented}.

The system consists of two different participants in the system exhibiting a trade-off relationship. In this case, we name them: prey and predator. The cycle dynamics start from the increasing population of prey. Then, the abundance of prey would increase the population of the predator. When the predator population increases, the population of prey decreases. Then due to a lack of food sources, the population of predators decreases which results in an increase in the population of prey.  
Their population are described in equation (1) and (2).


$$
x_{1t} = x_{1t}
$$

$$
\begin{gather}
    x_{1t} = x_{1t} [r_1 - a_{12}x_{2t}] \\
    x_{2t} = x_{2t} [-r_2 - a_{21}x_{1t}] 
\end{gather}
$$

In this mode, ${x}_{1t}$ and ${x}_{2t}$ are the population densities of prey and a predator at time t. The $\dot{x}_{1t}$ and $\dot{x}_{2t}$ are the change rate for the population densities, The $r_1$ and $r_2$ are the natural growth and death rate for the predator and prey, and $a_{ij}$ is the interaction parameter between different species.

We assume the growth rate and species interaction parameters are fixed in the system, and we take one step approximation of the above equation to the approximated next time step state of our system as shown in equation.

$$
\begin{gather}
    x_{1,t+1} = \Delta t \cdot \dot{x}_{1,t} + x_{1,t} + v_{1} \\
    x_{2,t+1} = \Delta t \cdot \dot{x}_{2,t} + x_{2,t} + v_{2} 
\end{gather}
$$
$$
\frac{1}{\Bigl(\sqrt{\phi \sqrt{5}}-\phi\Bigr) e^{\frac25 \pi}} \equiv 1+\frac{e^{-2\pi}} {1+\frac{e^{-4\pi}} {1+\frac{e^{-6\pi}} {1+\frac{e^{-8\pi}} {1+\cdots} } } }
$$
The $\Delta t$ here is the time step. Here, process noise $v_1$ and $v_2$ as zero mean and variance $\sigma_1$ and $\sigma_2$. Clearly, our system is not linear. Thus it can be represented in a matrix form for state transition. We rewrite the system in equation \ref{eq:fsys}.

$$
\begin{equation}
    \begin{bmatrix}
    x_{1,t+1}\\
    x_{2,t+1}
    \end{bmatrix}
    = f(\begin{bmatrix}
    x_{1,t}\\
    x_{2,t}
    \end{bmatrix}) + \begin{bmatrix}
    v_{1}\\
    v_{2}
    \end{bmatrix}
\end{equation}
$$

As for the measurement, we assume that we are able to measure the $x_{n,1}$ and $x_{n,2}$ directly with a certain variance. In this case, the measurement is linear. Then the measurement function can be written in equation \ref{eq:measurement}. The $u_1$ and $u_2$ denote measurement noise.

$$
\begin{equation}
    \begin{bmatrix}
    y_{1,t}\\
    y_{2,t}
    \end{bmatrix}
    = \begin{bmatrix}
    1 & 0\\
    0 & 1
    \end{bmatrix} \cdot \begin{bmatrix}
    x_{1,t}\\
    x_{2,t}
    \end{bmatrix} + \begin{bmatrix}
    u_{1}\\
    u_{2}
    \end{bmatrix}
\end{equation}
$$

## Linearized Kalman Filter


## Extended Kalman Filter


## Unscented Kalman Filter

# Simulation

## Simulation Setting

## Tracking Simulation

## Performance Benchmarking

# Conclusion 
The above simulation shows that all three algorithms can track the Lotka-Volterra system with some variance. Among all three algorithms, UKF tracking lowers RMSE but has the highest computation cost. The processing time from UKF is much larger than the processing time from EKF and LKF. The RMSE from EKF is slightly larger than UKF, but its processing time is much lower than UKF. In such a case, the EKF is a better choice among all these candidate algorithms for the Lotka-Volterra system since it offers an outstanding accuracy with low complexity.

The main reason for the processing time difference is the operation complexity of different algorithms. For the EKF, the most time-consuming operation is the Jacobian matrix formulation. As for the UKF, the most time-consuming operation is the Cholesky factorization during sampling for the covariance matrix. Since the system's order is only first order, the computation cost for Jacobian matrix formulation is low. Thus, in our case, the EKF's complexity is lower than the UKF. Also, in most literature, the UKF can offer better performance in terms of RMSE than the EKF. Depending on the system, the performance gain from UKF can be really. For some cases, trading for high accuracy with high complexity is more desirable. In such conditions, the UKF is more desirable than EKF.
Regarding stability, the UKF could fail during the iteration when performing Cholesky factorization. However, there are specific ways in literature to avoid this, such as iteratively factorization\cite{chandrasekar2008reduced}. However, this would downgrade the performance of UKF. Thus, in terms of stability, the EKF can outperform UKF.

The project can be extended to more complex system trackings with other non-linear filter algorithms for further work. For the Lotka-Volterra system, the tracking performance comparison is not that significant due to the complexity. To further benchmark, we can increase the system complexity. Also, there are other candidate algorithms for non-linear system tracking, such as the Particle filter. For further work, our goal is to implement other algorithms for comparison. 

# Further Readings