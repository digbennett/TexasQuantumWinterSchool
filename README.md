<table>
  <tr>
    <td align="center" valign="middle" style="vertical-align: middle; text-align: center; padding-right: 16px;">
      <img src="./TQWS_Background.jpg" alt="Texas Quantum Winter School" width="260" />
    </td>
    <td>
      <h1>Texas Quantum Winter School</h1>
      <p>
        This repository contains information and materials associated with the
        <strong>Texas Quantum Winter School</strong>, an intensive educational program focused on open
        quantum systems and quantum information in the chemical sciences.
      </p>
      <p><strong>Official website:</strong> <a href="https://www.texquantum.com">https://www.texquantum.com</a></p>
    </td>
  </tr>
</table>

<hr/>

## Program Description

The **Texas Quantum Winter School** is a week-long school designed to provide 
early-career researchers with a rigorous introduction to theoretical and computational 
approaches for open quantum systems, quantum dynamics, and quantum information concepts 
relevant to chemical physics.

The program emphasizes foundational theory, practical tools, and close interaction 
between participants and instructors.

- **Location:** Wimberley, Texas, USA  
- **Dates:** January 4–9, 2026  

The school is supported by the **National Science Foundation (NSF)**.

---

## Lecturers

Instruction is provided by researchers with expertise in open quantum systems, quantum 
information, and chemical physics.

* Prof. Ignacio Franco (U Rochester) - Open Quantum Systems
* Prof. Eric Bittner (U Houston) - Quantum Light
* Prof. Doran Raccah (U Texas, Austin) - Open Quantum Systems
* Prof. Ajay Kandada (Wake Forrest) - Quantum Light
* Prof. Gabriel Landi (U Rochester) - Quantum Information
* Prof. Kade Head-Marsden (U Minnesota) - Quantum Computing

---

## Talk Descriptions

* **Monday, January 5th**
  * Talk 1: Quantum Information Lecture 1 [Prof. Gabriel Landi]
    * Subtopic: Density Matrices
    * Description: This lecture introduces density matrices as a general framework for 
    describing quantum states, emphasizing how they encode all experimentally 
    accessible information and unify pure and mixed-state descriptions. It motivates 
    this perspective through measurement, statistical interpretation, and composite 
    quantum systems.
  * Talk 2: Quantum Information Lecture 2 [Prof. Gabriel Landi]
    * Subtopic: Quantum Channels
    * Description: This lecture develops a high-level perspective on quantum channels 
    as a way of thinking about how quantum states can evolve and be processed beyond 
    idealized, closed-system dynamics. The emphasis is on framing open evolution, 
    noise, and measurement as related manifestations of a common underlying structure, 
    rather than as separate phenomena.
  * Talk 3: Quantum Computing Lecture 1 [Prof. Kade Head-Marsden]
    * Subtopic: Spin-Spin Decoherence
    * Description: This lecture introduces Kraus operators to describe amplitude 
    damping, phase damping, and decoherence. Emphasis is given on the timescales of 
    each process and standard nomenclature.
  * Talk 4: Quantum Computing Lecture 2 [Prof. Kade Head-Marsden]
    * Subtopic: Adapting Central Spin Problem
    * Description: This lecture focuses on spin-spin relaxation in the central spin 
    model including methods to treat these dynamics such as the cluster correlation 
    expansion (CCE) and the analytical pair product approximation (APPA).


* **Tuesday, January 6th**
  * Talk 5: Quantum Information Lecture 3 [Prof. Gabriel Landi]
    * Subtopic: Lindblad Equations
    * Description: Derives the Lindblad (GKSL) master equation and clarifies complete 
    positivity, trace preservation, and the structure of dissipators.
  * Talk 6: Quantum Light Lecture 1 [Prof. Ajay Kandada]
    * Subtopic: Intro Non-Linear Optics
    * Description: This lecture introduces the quantum mechanical description of light 
    and photodetection, emphasizing how nonclassical states—such as squeezed and 
    entangled photons generated via nonlinear optical processes—enable measurements 
    beyond classical shot-noise limits.
  * Talk 7: Open Quantum Systems Lecture 1 [Prof. Ignacio Franco]
    * Subtopic: Quantum Master Equations (HEOM)
    * Description: Introduces the hierarchical equations of motion (HEOM) for 
    non-Markovian open quantum dynamics. Covers bath spectral densities, convergence 
    control, and applications to chemical excitations and charge transport.
  * Talk 8: Quantum Light Lecture 2 [Prof. Eric Bittner]
    * Subtopic: Stochastic Processes
    * Description: The discussion focuses on the physical meaning of Itō versus 
    Stratonovich calculus and on how measurement backaction and dissipation lead to 
    effective violations of PT symmetry, clarifying when and why symmetry arguments 
    break down in realistic quantum optical systems.


* **Wednesday, January 7th**
  * Talk 9: Quantum Information Lecture 4 [Prof. Gabriel Landi]
    * Subtopic: Correlation/Entropy/Entanglement
    * Description: Defines and relates classical and quantum correlations via von 
    Neumann entropy, mutual information, and entanglement measures.
  * Talk 10: Open Quantum Systems Lecture 2 [Prof. Doran Raccah]
    * Subtopic: Stochastic Methods in a Coherent State Basis
    * Description: Derives the non-Markovian stochastic Schrödinger equation using 
    Bargmann coherent-states. Demonstrates the localization which arises in the wide 
    open quantum system model.
  * Talk 11: Open Quantum Systems Lecture 3 [Prof. Ignacio Franco]
    * Subtopic: HEOM and pytenso
    * Description: This lecture consisted of a hands-on session using tensor-network 
    HEOM tools (pytenso) to accelerate simulations. Demonstrates model setup, 
    truncation strategies, and analysis of dynamics with provided notebooks.
  * Talk 12: Quantum Light Lecture 3 [Prof. Eric Bittner & Prof. Ajay Kandada]
    * Subtopic: Beyond Classical Observables
    * Description: Includes topics such as topology and spectral entanglement in 
    cavity-mediated photon scattering.


* **Thursday, January 8th**
  * Talk 13: Quantum Information Lecture 5 [Prof. Gabriel Landi]
    * Subtopic: Quantum Metrology
    * Description: Studied transport across a Fibonacci chain by constructing a master 
    equation in Mathematica using Melt.
  * Talk 14: Quantum Computing Lecture 3 [Prof. Kade Head-Marsden]
    * Subtopic: Introduction to Quantum Circuits
    * Description: Covers the circuit model, single- and two-qubit gates, measurement, 
    and basic circuit identities. Builds intuition for composing circuits and analyzing 
    noise and depth/width trade-offs.
  * Talk 15: Open Quantum Systems Lecture 4 [Prof. Doran Raccah]
    * Subtopic: HOPS
    * Description: Introduces the hierarchy of pure states (HOPS) method for 
    non-Markovian dynamics and its relation to HEOM. Details numerical construction, 
    convergence, and adaptivity, along with representative applications.
  * Talk 16: Quantum Computing Lecture 4 [Prof. Kade Head-Marsden]
    * Subtopic: Algorithms for Energy Estimation
    * Description: Surveys energy-estimation algorithms including VQE and 
    phase-estimation variants. Emphasizes ansatz design, error mitigation, and 
    applications to molecular electronic structure.

## Repository Contents
Top-level folders are organized by topic; each topic contains per‑lecture subfolders 
("Lecture N: Title") with slides, exercises, and notebooks/code. Shared files are 
stored at the repository root.

Structure at a glance:
```
TexasQuantumWinterSchool/
├── Topic/
│   └── Lecture N: Title/
│       ├── slides/
│       ├── exercises/
│       └── notebooks/ or code/
└── Root files: README.md, LICENSE, TQWS_Background.jpg
```

---

## Contact

For additional information, please visit the official website or contact:

**Email:** texasquantumwinterschool@gmail.com
