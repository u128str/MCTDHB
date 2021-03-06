
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
<meta name="generator" content="http://www.nongnu.org/elyxer/"/>
<meta name="create-date" content="2015-04-01"/>
<link rel="stylesheet" href="http://elyxer.nongnu.org/lyx.css" type="text/css" media="all"/>
<title>Exact quantum dynamics of bosons with finite-range time-dependent interactions of harmonic type</title>
</head>
<body>
<div id="globalWrapper">
<div class="Unindented">
The main goal of the MCTDHB-Lab project with Template.1 (default) is to reproduce some results from Table <a class="Reference" href="#tabHIMenergies">1↓</a> of the manuscript (see below) published in <a class="URL" href="http://link.aps.org/doi/10.1103/PhysRevA.86.063606">PRA 86, 063606 (2012)</a>:
</div>
<div class="Indented">
<div class="listing">
<pre class="listing">1)   Verify that: 
    all the parameters of the Hamiltonian are correct
     click on "Hamitonian"-&gt;"System"-&gt;"N"-button (similarly Trap V(x), strength of the interaction lambda, etc.)
</pre>
</div>
<img class="embedded" src="step_1.png" alt="figure step_1.png" style="max-width: 1022px; max-height: 594px;"/>
 
</div>
<div class="Indented">
<div class="listing">
<pre class="listing">2)   Verify that: 
     TDSE is solved for relaxation, i.e., you are looking for the ground state
     click on "Hamiltonian"-&gt;"TDSE" push on the "Relaxation"-button
</pre>
</div>
<img class="embedded" src="step_1b.png" alt="figure step_1b.png" style="max-width: 929px; max-height: 196px;"/>

</div>
<div class="Indented">
<div class="listing">
<pre class="listing">3)   Run MCTDHB (click on "Computation"-&gt;"Run MCTDHB"-&gt;"Run"-&gt;"Ok"-button), wait till it is finished 
     click on "Refresh"-button to see outcome of the current iteration
     click on "Proj.Size"-button in the left corner of the MCTDHB-Lab windows
     to see how long one has to wait untill the end of the current run
</pre>
</div>
<img class="embedded" src="step_2.png" alt="figure step_2.png" style="max-width: 1434px; max-height: 775px;"/>

</div>
<div class="Indented">
<div class="listing">
<pre class="listing">4)   Compare the obtained energy (from the last printed iteration) with the analytical results depicted in the Table
5)   Change level M of the MCTDHB(M) method:  "Hamiltonian"-&gt;"Morb and Psi at t=0"-&gt;"M"-button
</pre>
</div>
<img class="embedded" src="step_4.png" alt="figure step_4.png" style="max-width: 1234px; max-height: 224px;"/>

</div>
<div class="Indented">
<div class="listing">
<pre class="listing">5) Run MCTDHB again, compare the results
</pre>
</div>

</div>
<div class="Indented">
<div class="float">
<a class="Label" name="Table-1"> </a><div class="table">
<div class="center">
<table>
<tr>
<td align="center" valign="top">
M 
</td>
<td align="center" valign="top">
N=10 
</td>
<td align="center" valign="top">
N=100 
</td>
<td align="center" valign="top">
N=1000 
</td>

</tr>
<tr>
<td align="center" valign="top">
1 
</td>
<td align="center" valign="top">
7.0<i><span class="red"> 71067811865483</span></i>
</td>
<td align="center" valign="top">
70.<span class="red"> 71067811865483</span> 
</td>
<td align="center" valign="top">
707.<span class="red"> 1067811865483</span> 
</td>

</tr>
<tr>
<td align="center" valign="top">
2 
</td>
<td align="center" valign="top">
7.038<span class="red"> 769026303168</span>
</td>
<td align="center" valign="top">
70.6801<span class="red"> 6951747168</span> 
</td>
<td align="center" valign="top">
707.0764<span class="red"> 334257315</span> 
</td>

</tr>
<tr>
<td align="center" valign="top">
3 
</td>
<td align="center" valign="top">
7.0383<span class="red"> 50652406389</span>
</td>
<td align="center" valign="top">
70.680125<span class="red"> 41218675</span> 
</td>
<td align="center" valign="top">
707.07642898<span class="red"> 71865</span> 
</td>

</tr>
<tr>
<td align="center" valign="top">
4 
</td>
<td align="center" valign="top">
7.0383484<span class="red"> 24909910</span>
</td>
<td align="center" valign="top">
70.6801253917<span class="red"> 4549</span> 
</td>
<td align="center" valign="top">

</td>

</tr>
<tr>
<td align="center" valign="top">
5 
</td>
<td align="center" valign="top">
7.0383484153<span class="red"> 49058</span>
</td>
<td align="center" valign="top">
70.680125391737<span class="red"> 62</span> 
</td>
<td align="center" valign="top">

</td>

</tr>
<tr>
<td align="center" valign="top">
6 
</td>
<td align="center" valign="top">
7.038348415311<span class="red"> 494</span>
</td>
<td align="center" valign="top">

</td>
<td align="center" valign="top">

</td>

</tr>
<tr>
<td align="center" valign="top">
7 
</td>
<td align="center" valign="top">
7.03834841531101 <span class="red">8</span>
</td>
<td align="center" valign="top">

</td>
<td align="center" valign="top">

</td>

</tr>
<tr>
<td align="center" valign="top">
<span class="formula"><i>E</i><sub><span class="mathrm">exact</span></sub></span> 
</td>
<td align="center" valign="top">
7.038348415311011 
</td>
<td align="center" valign="top">
70.68012539173752 
</td>
<td align="center" valign="top">
707.0764289869851 
</td>

</tr>

</table>

</div>
<div class="caption">
Table 1 Ground state energies of the harmonic interaction Hamiltonian for the systems of N=10,100,1000 bosons. Exact analytical versus numerical MCTDHB(<span class="formula"><i>M</i></span>) results, <span class="formula"><i>M</i></span> is the number of self-consistent orbitals used. The interparticle interaction strengths have been chosen to keep <span class="formula">Λ = <i>K</i><sub>0</sub>(<i>N</i> − 1) = 0.5</span> constant. In this case all these systems have the same Gross-Pitaevskii solution, i.e., the same energy per particle. The one orbital MCTDHB(<span class="formula"><i>M</i> = 1</span>) theory is fully equivalent to the Gross-Pitaevskii mean-field. It is seen that converged results are obtained with less self-consistent orbitals when increasing the number of particles.
</div>
<a class="Label" name="tabHIMenergies"> </a> 
</div>

</div>

</div>
<h1 class="title">
Exact quantum dynamics of bosons with finite-range time-dependent interactions of harmonic type
</h1>
<h2 class="author">
Axel U. J. Lode<span class="formula"><sup>1</sup></span>, Kaspar Sakmann<span class="formula"><sup>1</sup></span>, Ofir E. Alon<span class="formula"><sup>2</sup></span>, Lorenz S. Cederbaum<span class="formula"><sup>1</sup></span>, and Alexej I. Streltsov<span class="formula"><sup>1</sup></span><span class="FootOuter"><span class="SupFootMarker"> [A] </span><span class="HoverFoot"><span class="SupFootMarker"> [A] </span>E-mail: Alexej.Streltsov@pci.uni-heidelberg.de</span></span> 
</h2>
<div class="Affiliation">
<span class="formula"><sup>1</sup></span> Theoretische Chemie, Physikalisch-Chemisches Institut, Universität Heidelberg,<br/>
 Im Neuenheimer Feld 229, D-69120 Heidelberg, Germany
</div>
<div class="Affiliation">
<span class="formula"><sup>2</sup></span> Department of Physics, University of Haifa at Oranim, Tivon 36006, Israel
</div>
<div class="abstract">
<p class="abstract-message">
Abstract
</p>
The exactly solvable quantum many-particle model with harmonic one- and two-particle interaction terms is extended to include time-dependency. We show that when the external trap potential and finite-range interparticle interaction have a time-dependency the exact solutions of the corresponding time-dependent many-boson Schrödinger equation are still available. We use these exact solutions to benchmark the recently developed multiconfigurational time-dependent Hartree method for bosons (MCTDHB) [Phys. Rev. Lett. <b>99</b>, 030402 (2007), Phys. Rev. A <b>77</b>, 033613 (2008)]. In particular, we benchmark the MCTDHB method for: (i) the ground state; (ii) the breathing many-body dynamics activated by a quench scenario where the interparticle interaction strength is suddenly turned on to a finite value; (iii) the non-equilibrium dynamic for driven scenarios where both the trap- and interparticle-interaction potentials are <i>time-dependent</i>. Excellent convergence of the ground state and dynamics is demonstrated. The great relevance of the self-consistency and time-adaptivity, which are the intrinsic features of the MCTDHB method, is demonstrated by contrasting the MCTDHB predictions and those obtained within the standard full configuration interaction method spanning the Fock space of the same size, but utilizing as one-particle basis set the fixed-shape eigenstates of the one-particle potential. Connections of the model’s results to ultra-cold Bose-Einstein condensed systems are addressed. 
</div>
<div class="PACS">
03.75.Kk, 05.30.Jp, 03.65.-w
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-1">1</a> Introduction
</h1>
<div class="Unindented">
Since the first realizations of Bose-Einstein condensates (BECs) <span class="bibcites">[<a class="bibliocite" name="cite-42" href="#biblio-42"><span class="bib-index">42</span></a>, <a class="bibliocite" name="cite-12" href="#biblio-12"><span class="bib-index">12</span></a>, <a class="bibliocite" name="cite-9" href="#biblio-9"><span class="bib-index">9</span></a>]</span> the experiments on this unique state of quantum systems have become more and more complex. Nowadays, the strength of the interparticle interactions, the trapping potential and the dimensionality of BECs are under experimental control <span class="bibcites">[<a class="bibliocite" name="cite-2" href="#biblio-2"><span class="bib-index">2</span></a>, <a class="bibliocite" name="cite-33" href="#biblio-33"><span class="bib-index">33</span></a>, <a class="bibliocite" name="cite-14" href="#biblio-14"><span class="bib-index">14</span></a>, <a class="bibliocite" name="cite-11" href="#biblio-11"><span class="bib-index">11</span></a>, <a class="bibliocite" name="cite-17" href="#biblio-17"><span class="bib-index">17</span></a>]</span>. This makes BECs a vivid and rich testing ground for a wide range of physical theories. Recent realizations of the dipolar BECs <span class="bibcites">[<a class="bibliocite" name="cite-3" href="#biblio-3"><span class="bib-index">3</span></a>, <a class="bibliocite" name="cite-21" href="#biblio-21"><span class="bib-index">21</span></a>, <a class="bibliocite" name="cite-26" href="#biblio-26"><span class="bib-index">26</span></a>]</span> open a new perspective in the development of the physics of ultra-cold atoms and molecules. A control on a new degree-of-freedom is achieved — the dipolar long-range part of the interparticle interaction can be now customized. This achievement can be considered as a first successful step towards a control on the overall shape of the interparticle interaction. It also stimulates the development of theoretical methods capable to solve the time-(in)dependent many-particle Schrödinger equation (TDSE) which governs the physics of the trapped ultra-cold systems with general interparticle interactions. The class of the many-body Hamiltonians which permits analytical solutions is quite small, so, generally, one has to rely on numerical many-body methods to solve the TDSE. The many-body methods in use have to be qualified to describe quantum many-body statics and dynamics. Benchmarking of these methods against exactly-solvable Hamiltonians is a necessary step for such a qualification.
</div>
<div class="Indented">
In this work we consider an exactly solvable many-body Hamiltonian, where both the one-body (trap) and two-body (interparticle interaction) potentials are of harmonic type, also known as the harmonic interaction model (HIM), see Refs. <span class="bibcites">[<a class="bibliocite" name="cite-25" href="#biblio-25"><span class="bib-index">25</span></a>, <a class="bibliocite" name="cite-23" href="#biblio-23"><span class="bib-index">23</span></a>]</span>. The exact solutions of the HIM problem are obtained by transformation of the Hamiltonian from the laboratory to the center of mass frame, where the Hamiltonian becomes separable. The price of the transformation is that an intuitive physical picture of &ldquo;real&rdquo; particles is lost and, instead, one operates with effective &ldquo;particles&rdquo; representing the transformed coordinates. One wants to have, first, a general many-body method for identical physical particles where each particle has its own &ldquo;real&rdquo; coordinate. Second, the method must be powerful enough to solve problems where such &ldquo;real&rdquo; coordinates are not favorable (suitable). Third, it should be capable to solve general problems where separations (transformations) of the variables are impossible, as it is the case in unharmonic and multiwell traps. In this work we want to test such a method — the recently developed multiconfigurational time-dependent Hartree for bosons (MCTDHB) method <span class="bibcites">[<a class="bibliocite" name="cite-35" href="#biblio-35"><span class="bib-index">35</span></a>, <a class="bibliocite" name="cite-29" href="#biblio-29"><span class="bib-index">29</span></a>]</span>, which can treat the dynamics, i.e., the TDSE for trapped bosons with general interparticle interactions. We would like to carry the examples to an extreme case separable in suitable coordinates; in this way we can unambiguously test the performance of the MCTDHB against analytical/exact solutions.
</div>
<div class="Indented">
While there are several known many-body models with time-independent Hamiltonians which have analytical solutions, the exactly solvable time-dependent many-body Hamiltonians are even less abundant. Why one needs to study them? Apart from an exploration of novel dynamical physical phenomena, there is a practical reason for it. In typical experiments with ultra-cold systems the manipulations of the trapping potentials, as well as altering of the magnetic field used in the Feshbach resonance technique(s) <span class="bibcites">[<a class="bibliocite" name="cite-11" href="#biblio-11"><span class="bib-index">11</span></a>]</span> to manipulate the interparticle interaction, are time-dependent procedures. So, there is a need for a proven theoretical method capable to solve time-dependent Hamiltonians where both the trap and interparticle interaction potentials are time-dependent. However, the TDSE with general time-dependent Hamiltonians can be solved only numerically, hence it is very difficult to verify and quantify the region of applicability and quality of the numerical solutions obtained. Convincing comparisons/benchmarks against exact results are of great relevance. In this work we show how to extend the exactly solvable quantum many-particle HIM problem to include time-dependency, and use it to benchmark the MCTDHB method.
</div>
<div class="Indented">
It is worthwhile to mention that some physical phenomena and properties of the many-body solutions of the HIM problem are &ldquo;universal&rdquo;, i.e., transferable to systems with other interparticle interactions, e.g., contact interaction. For example, small displacements of the density out of the center of an harmonic trap result in so called &ldquo;dipole oscillations&rdquo; with the trap frequency which are independent of the interparticle interaction. Another example is a quench of the interparticle interaction in an harmonically trapped system — it activates only &ldquo;breathing&ldquo; excitations which preserve the symmetry of the trap. In the present work we discover a novel time-dependent phenomenon in the extended HIM, and discuss its &ldquo;universality&rdquo; for the harmonically trapped systems with general interparticle interactions.
</div>
<div class="Indented">
The structure of the paper is as follows. In Sec. <a class="Reference" href="#HIMdfns">2↓</a> we introduce the harmonic interaction model and discuss the aftermaths and implications appearing due to the transformation of the coordinates from the center of mass frame, where the exact solutions are analytically known, to the laboratory frame where we want to solve the problem numerically. The MCTDHB method is briefly reviewed in Sec. <a class="Reference" href="#MCTDHB_th">3↓</a>. Sec. <a class="Reference" href="#GS">4↓</a> provides detailed benchmarks and comparisons of the exact and numerical results for the ground state of the HIM obtained within the framework of the MCTDHB and the standard full configuration interaction (FCI) methods, the latter is also known as the exact diagonalization (ED) technique. In Sec. <a class="Reference" href="#HIMquench">5↓</a> we benchmark our numerical tools to describe the breathing many-body dynamics activated by a quench scenario where the interparticle interaction strength is suddenly turned on from zero to a finite value. Sec. <a class="Reference" href="#TDHIM">6↓</a> shows how to extend the exactly solvable quantum many-particle model with harmonic one- and two-particle interaction terms to include time-dependency. Here we also demonstrate the applicability of the MCTDHB method to describe numerically-exact many-boson dynamics for complicated scenarios where both the external trap and interparticle interaction potentials are time-dependent. Sec. <a class="Reference" href="#sum">7↓</a> summarizes our results and outlooks the novel predictions obtained for the HIM problem to ultra-cold atomic systems with contact interactions.
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-2">2</a> The Harmonic Interaction Model (HIM)
</h1>
<div class="Unindented">
<a class="Label" name="HIMdfns"> </a> 
</div>
<h2 class="Subsection">
<a class="toc" name="toc-Subsection-2.1">2.1</a> Basic Definitions
</h2>
<div class="Unindented">
<a class="Label" name="HIMdf_A"> </a>
</div>
<div class="Indented">
Our starting point is the harmonic interaction model (HIM), see, e.g., Refs. <span class="bibcites">[<a class="bibliocite" name="cite-25" href="#biblio-25"><span class="bib-index">25</span></a>, <a class="bibliocite" name="cite-23" href="#biblio-23"><span class="bib-index">23</span></a>]</span>. The Hamiltonian of the HIM is readily obtained in the laboratory frame of reference by setting the interparticle interaction potential <span class="formula"><i>Ŵ</i></span> and the one-body potential <span class="formula"><i>V̂</i></span> in the many-body Hamiltonian in dimensionless units, <div class="formula">
<a class="eqnumber" name="HAM">(1) </a><i>Ĥ</i> = <span class="limits"><sup class="limit"><i>N</i></sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> = 1</sub></span><span class="symbol">(</span><i>T̂</i>(<i>r⃗</i><sub><i>i</i></sub>) + <i>V̂</i>(<i>r⃗</i><sub><i>i</i></sub>)<span class="symbol">)</span> + <span class="limits"><sup class="limit"><i>N</i></sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> &lt; <i>j</i></sub></span><i>Ŵ</i>(<i>r⃗</i><sub><i>i</i></sub>, <i>r⃗</i><sub><i>j</i></sub>), 
</div>
to be harmonic: <div class="formula">
<a class="eqnumber" name="HIMint">(2) </a><i>Ŵ</i>(<i>r⃗</i><sub><i>i</i></sub>, <i>r⃗</i><sub><i>j</i></sub>) = <i>K</i><sub>0</sub><span class="symbol">(</span><i>r⃗</i><sub><i>i</i></sub> − <i>r⃗</i><sub><i>j</i></sub><span class="symbol">)</span><sup>2</sup>,  <i>V</i>(<i>r⃗</i>) = <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i><sup>2</sup><i>r⃗</i><sup>2</sup>.
</div>
Here, <span class="formula"><i>K</i><sub>0</sub></span> accounts for the strength of the two-body interaction and <span class="formula"><i>T̂</i>(<i>r⃗</i>) =  − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>r⃗</i></sub></span></span> is the kinetic energy operator. A positive value of <span class="formula"><i>K</i><sub>0</sub></span> corresponds to an attraction while a negative value means repulsion. In the case of a parabolic trapping potential, it is easy to see that the system becomes unbound when the value of <span class="formula"><i>K</i><sub>0</sub></span> is negative and big enough for the two-body repulsion to overcome the one-body harmonic trapping, i.e. <span class="formula"><i>K</i><sub>0</sub> &lt;  − <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i><sup>2</sup></span><span class="ignored">)/(</span><span class="denominator">2<i>N</i></span><span class="ignored">)</span></span></span>.
</div>
<div class="Indented">
Following Cohen and Lee in Ref. <span class="bibcites">[<a class="bibliocite" name="cite-25" href="#biblio-25"><span class="bib-index">25</span></a>]</span>, the Hamiltonian, Eqs. (<a class="Reference" href="#HAM">1↑</a>,<a class="Reference" href="#HIMint">2↑</a>), can be separated into <span class="formula"><i>N</i></span> independent harmonic oscillators by the following coordinate transformations: <div class="formula">
<a class="eqnumber" name="eq-3">(3) </a><span class="environment"><span class="arrayrow">
<span class="arraycell align-r">
<i>q⃗</i><sub><i>j</i></sub>
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>j</i>(<i>j</i> + 1)</span><span class="ignored">)</span></span></span><span class="ignored">)</span></span><span class="limits"><sup class="limit"><i>j</i></sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> = 1</sub></span>(<i>r⃗</i><sub><i>j</i> + 1</sub> − <i>r⃗</i><sub><i>i</i></sub>),   <i>j</i> = 1, ..., <i>N</i> − 1, 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
<i>q⃗</i><sub><i>N</i></sub>
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>N</i></span><span class="ignored">)</span></span></span><span class="ignored">)</span></span><span class="limits"><sup class="limit"><i>N</i></sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> = 1</sub></span><i>r⃗</i><sub><i>i</i></sub>.
</span>

</span>
</span>
</div>
The transformed Hamiltonian in the center of mass frame reads: <div class="formula">
<a class="eqnumber" name="HIMTD">(4) </a><span class="environment"><span class="arrayrow">
<span class="arraycell align-r">
<i>Ĥ</i>
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<i>Ĥ</i><sub><i>rel</i></sub> + <i>Ĥ</i><sub><i>CM</i></sub>, 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
<i>Ĥ</i><sub><i>rel</i></sub>
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="limits"><sup class="limit"><i>N</i> − 1</sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> = 1</sub></span>( − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>q⃗</i><sub><i>i</i></sub></sub></span> + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>δ</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>N</i></sub></span><i>q⃗</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>i</i></sub></span>),   
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
  <i>Ĥ</i><sub><i>CM</i></sub>
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
 − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>q⃗</i><sub><i>N</i></sub></sub></span> + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i><sup>2</sup><i>q⃗</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>N</i></sub></span>.
</span>

</span>
</span>
</div>
Here, <span class="formula"><i>δ</i><sub><i>N</i></sub> = <span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><sup>2</sup> + 2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span></span> is the trapping frequency of the <span class="formula"><i>N</i> − 1</span> harmonic oscillators originating from the set of relative coordinates; and <span class="formula">1</span> harmonic oscillator with the frequency <span class="formula"><i>ω</i></span> representing a center of mass coordinate. This separability of the HIM Hamiltonian into the center of mass and relative coordinates allows the following visualization: the overall HIM system can be pictured as a medium formed by <span class="formula"><i>N</i> − 1</span> identical, noninteracting particles associated with relative coordinates <span class="formula"><i>q</i><sub><i>k</i></sub>, <i>k</i> = 1, ..., <i>N</i> − 1</span>, moving in an effective harmonic trap with a time-independent frequency <span class="formula"><i>δ</i><sub><i>N</i></sub></span>, and an independent effective particle with coordinate <span class="formula"><i>q</i><sub><i>N</i></sub></span>, representing the system’s center of mass, trapped in the original time-independent harmonic potential with frequency <span class="formula"><i>ω</i></span>.
</div>
<div class="Indented">
The general solution of the HIM problem in its separable form, Eq. (<a class="Reference" href="#HIMTD">4↑</a>), is a product of <span class="formula"><i>N</i></span> generally different harmonic oscillator wavefunctions, and the total energy is the sum of the corresponding oscillator’s energies. The exact energy <span class="formula"><i>E</i><sub><i>exact</i></sub></span> of the ground state takes on a very simple form, see, e.g., Refs. <span class="bibcites">[<a class="bibliocite" name="cite-25" href="#biblio-25"><span class="bib-index">25</span></a>, <a class="bibliocite" name="cite-23" href="#biblio-23"><span class="bib-index">23</span></a>]</span>: <div class="formula">
<a class="eqnumber" name="himEexact">(5) </a><span class="environment"><span class="arrayrow">
<span class="arraycell align-r">
<i>E</i><sub><i>exact</i></sub>
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="fraction"><span class="ignored">(</span><span class="numerator"><i>D</i></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>N</i> − 1)<i>δ</i><sub><i>N</i></sub> + <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>D</i></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i>.
</span>

</span>
</span>
</div>
Here, <span class="formula"><i>D</i></span> is the dimensionality of the HIM system. We note that the HIM problem is an example of a many-body system with finite-range interparticle interactions which permits analytical solution in any dimension, i.e., in 1D and in higher dimensions. This is another attractive feature of the HIM model relevant for benchmarking numerical methods for the many-particle TDSE.
</div>
<h2 class="Subsection">
<a class="toc" name="toc-Subsection-2.2">2.2</a> Representing the HIM with Basis Functions
</h2>
<div class="Unindented">
Let us consider the ground state of the (one-dimensional) HIM problem for <span class="formula"><i>N</i> = 2</span> bosons. We contrast the solution written in the center of mass frame, i.e., in the <span class="formula"><i>q</i><sub>1</sub>, <i>q</i><sub>2</sub></span> coordinates with that in the laboratory frame, i.e., in <span class="formula"><i>x</i><sub>1</sub>, <i>x</i><sub>2</sub></span> coordinates. The corresponding transformations of the coordinates are given in Eq. (<a class="Reference" href="#coordtrans">↓</a>) with <span class="formula"><i>q⃗</i><sub><i>j</i></sub> = <i>q</i><sub><i>j</i></sub></span> and <span class="formula"><i>r⃗</i><sub><i>j</i></sub> = <i>x</i><sub><i>j</i></sub></span>. It is convenient to denote the <span class="formula"><i>n</i></span>-th harmonic oscillator (HO) function as <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>n</i></sub></span>(<i>X</i>, Ω)</span>, here <span class="formula">Ω</span> is harmonic oscillator frequency and <span class="formula"><i>X</i></span> a general variable. The ground state solution of the two-particle HIM problem reads: <div class="formula">
<a class="eqnumber" name="HIMGS2bosons">(6) </a><span class="environment"><span class="arrayrow">
<span class="arraycell align-r">
Ψ(<i>q</i><sub>1</sub>, <i>q</i><sub>2</sub>)
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>q</i><sub>1</sub>, <i>δ</i><sub>2</sub>)<i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>q</i><sub>2</sub>, <i>ω</i>) = <span class="scriptfont">N</span><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>δ</i><sub>2</sub><i>q</i><span class="scripts"><sup class="script">2</sup><sub class="script">1</sub></span></sup><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i><i>q</i><span class="scripts"><sup class="script">2</sup><sub class="script">2</sub></span></sup> ≡ 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
Ψ(<i>x</i><sub>1</sub>, <i>x</i><sub>2</sub>)
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="scriptfont"> N</span><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>δ</i><sub>2</sub>(<span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root">2</span><span class="ignored">)</span></span></span><span class="ignored">)</span></span>(<i>x</i><sub>2</sub> − <i>x</i><sub>1</sub>))<sup>2</sup></sup><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i>(<span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root">2</span><span class="ignored">)</span></span></span><span class="ignored">)</span></span>(<i>x</i><sub>2</sub> + <i>x</i><sub>1</sub>))<sup>2</sup></sup>
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">

</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="scriptfont"> N</span><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>δ</i><sub>2</sub> + <i>ω</i></span><span class="ignored">)/(</span><span class="denominator">4</span><span class="ignored">)</span></span><i>x</i><span class="scripts"><sup class="script">2</sup><sub class="script">2</sub></span></sup><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>δ</i><sub>2</sub> + <i>ω</i></span><span class="ignored">)/(</span><span class="denominator">4</span><span class="ignored">)</span></span><i>x</i><span class="scripts"><sup class="script">2</sup><sub class="script">1</sub></span></sup><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>ω</i> − <i>δ</i><sub>2</sub>)<i>x</i><sub>2</sub><i>x</i><sub>1</sub></sup>
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">

</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="scriptfont"> N</span><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>x</i><sub>1</sub>, <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>)<i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>x</i><sub>2</sub>, <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>)<i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>ω</i> − <i>δ</i><sub>2</sub>)<i>x</i><sub>2</sub><i>x</i><sub>1</sub></sup>
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">

</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="scriptfont"> N</span><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>x</i><sub>1</sub>, <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>)<i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>x</i><sub>2</sub>, <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>)<span class="array"><span class="arrayrow"><span class="bracket align-left">⎛</span></span><span class="arrayrow"><span class="bracket align-left">⎝</span></span></span>1 − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>ω</i> − <i>δ</i><sub>2</sub>)<i>x</i><sub>2</sub><i>x</i><sub>1</sub> + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">8</span><span class="ignored">)</span></span>(<i>ω</i> − <i>δ</i><sub>2</sub>)<sup>2</sup><i>x</i><span class="scripts"><sup class="script">2</sup><sub class="script">2</sub></span><i>x</i><span class="scripts"><sup class="script">2</sup><sub class="script">1</sub></span> − …<span class="array"><span class="arrayrow"><span class="bracket align-right">⎞</span></span><span class="arrayrow"><span class="bracket align-right">⎠</span></span></span>
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">

</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="limits"><sup class="limit">∞</sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> ≥ <i>j</i> = 0</sub></span><i>a</i><sub><i>ij</i></sub><i><span class="scriptfont">S</span>̂</i><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>i</i></sub></span>(<i>x</i><sub>1</sub>, <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>)<i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>j</i></sub></span>(<i>x</i><sub>2</sub>, <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>).
</span>

</span>
</span>
</div>
Here we use the Taylor expansion for the cross-term <span class="formula"><i>e</i><sup> − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>ω</i> − <i>δ</i><sub>2</sub>)<i>x</i><sub>2</sub><i>x</i><sub>1</sub></sup></span>. <span class="formula"><i><span class="scriptfont">S</span>̂</i></span> is the symmetrization operator and <span class="formula"><i>a</i><sub><i>ij</i></sub></span> — known reexpansion coefficients. The close inspection of the above transformation shows that a single Hartree product of two HO wavefunctions written in the center of mass frame is represented by an infinite sum of different Hartree products in the laboratory frame even if one uses the &ldquo;dressed&rdquo; frequencies <span class="formula"><span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span></span>. For bosonic systems these Hartree products have to be properly symmetrized to take into account permutational symmetry of the total wavefunction. Hence, a numerical solution of the HIM Hamiltonian in the laboratory frame is a very involved problem — the numerical convergence depends on how fast and efficient this sum is spanned.
</div>
<div class="Indented">
For numerical treatments the infinite sums, like in (<a class="Reference" href="#HIMGS2bosons">6↑</a>), must be truncated. The number of the terms <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub></span> considered defines the size of the Fock space spanned, i.e., the size of the respective secular matrix to be diagonalized in order to find the respective eigenvalues and eigenstates. For a general <span class="formula"><i>N</i></span>-boson system this size is <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub> = <span class="bigsymbol">(</span><span class="binom"><span class="binomstack"><i>N</i> + <i>M</i> − 1</span><span class="binomstack"><i>N</i></span></span><span class="bigsymbol">)</span></span>, where <span class="formula"><i>M</i></span> is the number of one-particle functions (orbitals) used to build the symmetrized Hartree products. The two-boson problem can be diagonalized by taking a lot of basis functions, while already for a ten-boson problem with <span class="formula"><i>M</i> = 16</span>, <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub> = 3, 268, 760</span>, meaning that the diagonalization of the respective secular matrix is not a simple task. Due to this binomial dependency of the size of the spanned Fock subspace, the bosonic system with large number of particles can be tackled only with quite a few orbitals, for example for <span class="formula"><i>N</i> = 1000</span> and <span class="formula"><i>M</i> = 3</span> the size of the secular problem is <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub> = 501, 501</span>, while already for <span class="formula"><i>M</i> = 4</span> it is <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub> = 167, 668, 501</span>. One of the main goals of the present work is to verify that the choice of the basis functions (orbitals) used to build up the permanents (symmetrized Hartree products) has enormous impact on the numerical convergence of many-body problems.
</div>
<div class="Indented">
The above considered analysis of the interplay between the exact solution written in the center of mass and in the laboratory frames of references is of course applicable for an arbitrary number of particles. The Hartree products appearing in the solutions of the corresponding HIM problems written in the laboratory frame are build up of HO basis functions with exponents <span class="formula"> − <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub><i>N</i></sub></span><span class="ignored">)/(</span><span class="denominator">4</span><span class="ignored">)</span></span><i>x</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>j</i></sub></span></span>, i.e., they depend via <span class="formula"><i>δ</i><sub><i>N</i></sub></span> on the number of particles <span class="formula"><i>N</i></span> and on the interparticle interaction strength <span class="formula"><i>K</i><sub>0</sub></span>. Hence, to solve the HIM problem for different parameters <span class="formula"><i>N</i></span> and <span class="formula"><i>K</i><sub>0</sub></span> in the laboratory frame it is advantageous to use one-particle basis functions with different shapes (exponents). It is natural to ask the following questions: What to do in a general case, when the analytic solution of the problem is not available, i.e., which basis set to use? And how to find the best &ldquo;optimal&rdquo; self-consistent orbitals? The simplest answer is to use the one-particle functions (bare orbitals) of the studied system when the two-body interparticle interactions are switched off. These fixed-shape basis functions are obtained as solutions of the one-particle problem <span class="formula"><i>ĥ</i><i>ψ</i><sub><i>i</i></sub> = <i>e</i><sub><i>i</i></sub><i>ψ</i><sub><i>i</i></sub></span> with <span class="formula"><i>ĥ</i> = <i>T̂</i>(<i>r⃗</i>) + <i>V̂</i>(<i>r⃗</i>)</span>. In the studied HIM model it means to use the HO basis <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>n</i></sub></span>(<i>x</i>, Ω)</span> with trapping frequency <span class="formula">Ω = <i>ω</i></span>. The answer to the second question is also known — to use the recently developed MCTDHB method <span class="bibcites">[<a class="bibliocite" name="cite-35" href="#biblio-35"><span class="bib-index">35</span></a>, <a class="bibliocite" name="cite-29" href="#biblio-29"><span class="bib-index">29</span></a>]</span>, which utilizes the Dirac-Frenkel variational principle to determine the optimal shapes of the orbitals for time-dependent problems. In this work we will examine and contrast the performance of two many-body methods to attack the time-dependent and time-independent HIM problems — the full configuration interaction (exact diagonalization) method which utilizes as a basis set the solutions of the one-particle in the harmonic trap, and the self-consistent MCTDHB method.
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-3">3</a> The MCTDHB Method
</h1>
<div class="Unindented">
<a class="Label" name="MCTDHB_th"> </a>
</div>
<div class="Indented">
Let us briefly describe the MCTDHB theory, for complete derivations and some recent applications see the literature <span class="bibcites">[<a class="bibliocite" name="cite-35" href="#biblio-35"><span class="bib-index">35</span></a>, <a class="bibliocite" name="cite-29" href="#biblio-29"><span class="bib-index">29</span></a>, <a class="bibliocite" name="cite-4" href="#biblio-4"><span class="bib-index">4</span></a>, <a class="bibliocite" name="cite-5" href="#biblio-5"><span class="bib-index">5</span></a>, <a class="bibliocite" name="cite-30" href="#biblio-30"><span class="bib-index">30</span></a>, <a class="bibliocite" name="cite-32" href="#biblio-32"><span class="bib-index">32</span></a>, <a class="bibliocite" name="cite-30" href="#biblio-30">30</a>, <a class="bibliocite" name="cite-36" href="#biblio-36"><span class="bib-index">36</span></a>, <a class="bibliocite" name="cite-16" href="#biblio-16"><span class="bib-index">16</span></a>, <a class="bibliocite" name="cite-33" href="#biblio-33">33</a>, <a class="bibliocite" name="cite-15" href="#biblio-15"><span class="bib-index">15</span></a>, <a class="bibliocite" name="cite-8" href="#biblio-8"><span class="bib-index">8</span></a>]</span>. The MCTDHB method has been developed to solve the time-(in)dependent many-boson Schrödinger equation. It relies on a multiconfigurational ansatz for the wavefunction, i.e., <span class="formula">∣Ψ(<i>t</i>)⟩ = <span class="limits"><span class="limit">∑</span></span><sub><i>n⃗</i></sub><i>C</i><sub><i>n⃗</i></sub>(<i>t</i>)∣<i>n⃗</i>;<i>t</i>⟩</span>. The unknown <span class="formula"><i>C</i><sub><i>n⃗</i></sub>(<i>t</i>)</span> are called expansion coefficients. The permanents <span class="formula">∣<i>n⃗</i>;<i>t</i>⟩</span> are build as symmetrized Hartree products of <span class="formula"><i>N</i></span> unknown orthogonal one-particle functions. For <span class="formula"><i>M</i></span> orbitals the number of these permanents is equal to the number all possible permutations of <span class="formula"><i>N</i></span> particles over <span class="formula"><i>M</i></span> orbitals, namely <span class="formula"><span class="bigsymbol">(</span><span class="binom"><span class="binomstack"><i>N</i> + <i>M</i> − 1</span><span class="binomstack"><i>N</i></span></span><span class="bigsymbol">)</span></span>. It is noteworthy that <i>both</i> the coefficients <span class="formula"><i>C</i><sub><i>n⃗</i></sub>(<i>t</i>)</span> and the one-partice functions used to build the permanents are time-dependent, variationally optimized quantities, which are determined by solving the corresponding MCTDHB equations of motion <span class="bibcites">[<a class="bibliocite" name="cite-35" href="#biblio-35"><span class="bib-index">35</span></a>, <a class="bibliocite" name="cite-29" href="#biblio-29"><span class="bib-index">29</span></a>]</span>. These equations depend on the parameters of the Hamiltonian, on the number of particles as well as on the number of the one-particle basis functions used. For different evolution times the optimal orbitals have different shapes — this feature is called time-adaptivity.
</div>
<div class="Indented">
Within the MCTDHB method the time-independent variational solutions are obtained by propagating the MCTDHB equations in imaginary time, which is equivalent to solving the stationary problem variationally, as developed in the multiconfigurational Hartree method (MCHB) for bosons, see Ref. <span class="bibcites">[<a class="bibliocite" name="cite-34" href="#biblio-34"><span class="bib-index">34</span></a>]</span>. Hence, the static solutions we give here qualify as test suits for how the standard (time-independent) variational principle is handled numerically by MCTDHB. From now on we call the time-dependent variational MCTDHB solution &ldquo;time-adaptive&rdquo; to distinguish it from the &ldquo;self-consistent&rdquo; static, i.e., time-independent MCTDHB solution. If the one-particle functions used are not allowed for optimization, the MCTDHB method boils down to the standard full configuration interaction (exact diagonalization) method. Thus, one can consider the MCTDHB method as an exact diagonalization method with time-adaptive (self-consistent) orbitals. For a given number of orbitals the dimension of the secular problem involved for the FCI(ED) and MCTDHB computations is the same, <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub> = <span class="bigsymbol">(</span><span class="binom"><span class="binomstack"><i>N</i> + <i>M</i> − 1</span><span class="binomstack"><i>N</i></span></span><span class="bigsymbol">)</span></span>. If only one self-consistent orbital is considered the MCTDHB theory boils down to the famous Gross-Pitaevskii (GP) mean-field theory widely and often successfully used to describe static and dynamics of condensed bosonic systems <span class="bibcites">[<a class="bibliocite" name="cite-13" href="#biblio-13"><span class="bib-index">13</span></a>, <a class="bibliocite" name="cite-6" href="#biblio-6"><span class="bib-index">6</span></a>, <a class="bibliocite" name="cite-10" href="#biblio-10"><span class="bib-index">10</span></a>]</span>.
</div>
<div class="Indented">
At this point it is very important to stress that the MCTDHB and standard FCI methods used in this work operate with general Hamiltonians in the laboratory frame of reference, so, the separability of the HIM model is not taken into account. But this formulation of the methods allows one to attack general, i.e., inseparable problems as has been done in, e.g., Refs. <span class="bibcites">[<a class="bibliocite" name="cite-30" href="#biblio-30"><span class="bib-index">30</span></a>, <a class="bibliocite" name="cite-32" href="#biblio-32"><span class="bib-index">32</span></a>, <a class="bibliocite" name="cite-36" href="#biblio-36"><span class="bib-index">36</span></a>, <a class="bibliocite" name="cite-16" href="#biblio-16"><span class="bib-index">16</span></a>, <a class="bibliocite" name="cite-46" href="#biblio-46">46</a>, <a class="bibliocite" name="cite-15" href="#biblio-15"><span class="bib-index">15</span></a>, <a class="bibliocite" name="cite-8" href="#biblio-8"><span class="bib-index">8</span></a>]</span>.
</div>
<div class="Indented">
The MCTDHB equations of motion are solved numerically efficiently with the MCTDHB program package <span class="bibcites">[<a class="bibliocite" name="cite-4" href="#biblio-4"><span class="bib-index">4</span></a>]</span>. The current study relies on the propagation of the orbitals’ equations of motion with a shared-memory parallelized implementation of the Adams-Bashforth Moulton predictor corrector integrator <span class="bibcites">[<a class="bibliocite" name="cite-41" href="#biblio-41"><span class="bib-index">41</span></a>]</span> and the coefficients’ equations of motion with a hybridly OpenMP-MPI parallelized short iterative Lanczos algorithm <span class="bibcites">[<a class="bibliocite" name="cite-37" href="#biblio-37"><span class="bib-index">37</span></a>]</span>. As primitive basis functions representing the self-consistent (time-dependent) orbitals we use either the HO discrete variable representation <span class="bibcites">[<a class="bibliocite" name="cite-44" href="#biblio-44"><span class="bib-index">44</span></a>]</span> or the fast Fourier transform collocation method utilizing hybrid OpenMP-MPI parallelization, see <span class="bibcites">[<a class="bibliocite" name="cite-4" href="#biblio-4"><span class="bib-index">4</span></a>]</span>.
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-4">4</a> Ground state of the HIM: MCTDHB and FCI vs. exact solution
</h1>
<div class="Unindented">
<a class="Label" name="GS"> </a>
</div>
<div class="Indented">
We begin by benchmarking the MCTDHB and FCI methods against the ground state of the one-dimensional HIM. We consider systems of <span class="formula"><i>N</i> = 2, 10, 50, 100, 1000</span> bosons trapped in the parabolic trap potential <span class="formula"><i>V</i>(<i>x</i>) = <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>x</i><sup>2</sup></span> with the inter-boson interaction strengths selected to keep <span class="formula">Λ = <i>K</i><sub>0</sub>(<i>N</i> − 1) = 0.5</span> constant. Such a choice of the interaction strengths implies that all these systems have the same GP solution, i.e., are equivalent at the mean-field level of description. To find the properties of convergence of the MCTDHB and FCI methods towards the exact solution of the HIM, it is instructive to successively increase the number of orbitals, <span class="formula"><i>M</i></span>, used in the computation. In Fig. <a class="Reference" href="#HIM_MCTDHB">1↓</a> we plot the relative difference between the ground state eigenenergy and the corresponding exact energy <span class="formula">(<i>E</i><sub><span class="mathrm">MB</span></sub> − <i>E</i><sub><span class="mathrm">exact</span></sub>) ⁄ <i>E</i><sub><span class="mathrm">exact</span></sub></span> as a function of number of orbitals <span class="formula"><i>M</i></span> used. The self-consistent many-body (MB) MCTDHB(<span class="formula"><i>M</i></span>) results are plotted by open symbols, the corresponding fixed-orbital FCI<span class="formula">(<i>M</i>)</span> results are depicted by filled symbols.
</div>
<div class="Indented">
The key observation seen in Fig. <a class="Reference" href="#HIM_MCTDHB">1↓</a> is that the numerical results converge towards the exact ones with increasing number of the orbitals used. The performance of the self-consistent MCTDHB method, however, by far exceeds that of the fixed-orbital FCI. Note the logarithmic scale and number of decades spanned! The proper choice of the many-body basis set is very crucial — within the same size of the Fock subspace spanned (dimension of the secular matrix <span class="formula"><i>N</i><sub><span class="mathrm">conf</span></sub></span>) one can get an improvement of about six to eight orders of magnitude! The results prove that the exact solutions of the HIM can be obtained numerically using the MCTDHB method with just a few self-consistent orbitals.
</div>
<div class="Indented">
Another striking feature of the MCTDHB method is its performance for different particle numbers. The convergence is faster for larger particle numbers at fixed <span class="formula">Λ = <i>K</i><sub>0</sub>(<i>N</i> − 1)</span>. This is anticipated, because at the mean-field GP<span class="formula"> ≡ </span>MCTDHB(1) level the considered systems of bosons are equivalent and the GP solution of the static HIM problem tends to an exact one in the thermodynamic (<span class="formula"><i>N</i> → ∞</span>) limit. Nevertheless, for large, but finite number of bosons the MCTDHB significantly improves the GP description. For example, for N=1000 the relative differences to the exact energy obtained within the GP and MCTDHB(3) are <span class="formula"> ~ 10<sup> − 3</sup></span> and <span class="formula"> ~ 10<sup> − 11</sup></span> percents respectively, see Fig. <a class="Reference" href="#HIM_MCTDHB">1↓</a>. So, the self-consistency becomes more and more relevant for systems made of larger particle numbers. In contrast to the MCTDHB, the fixed-orbital FCI method does not show such a tendency. The low performance of the full configuration interaction method utilizing the bare HO orbital basis set is evident from the above done two-boson analysis in Eq. (<a class="Reference" href="#HIMGS2bosons">6↑</a>). Instead of the HO eigenfunctions of trap potential <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>n</i></sub></span>(<i>x</i>, Ω = <i>ω</i>)</span> one has to use the HO basis functions with &ldquo;modified&rdquo; frequency <span class="formula">Ω = <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub><i>N</i></sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span></span>. However, in the general case, when an analytical solution is unavailable the only strict way to find the &ldquo;proper&rdquo; basis set is to solve the MCTDHB(<span class="formula"><i>M</i></span>) equations, which determine variationally the optimal one-particle functions, see Ref. <span class="bibcites">[<a class="bibliocite" name="cite-34" href="#biblio-34"><span class="bib-index">34</span></a>]</span>.
</div>
<div class="Indented">
To highlight the convergence of the MCTDHB(<span class="formula"><i>M</i></span>) method with the number of orbitals <span class="formula"><i>M</i></span> used we present in table <a class="Reference" href="#tabHIMenergies">1↑</a> the total ground state energies of the above considered systems of N=10,100,1000 bosons with <span class="formula">Λ = <i>K</i><sub>0</sub>(<i>N</i> − 1) = 0.5</span>. The exact ground state energies are from Refs. <span class="bibcites">[<a class="bibliocite" name="cite-25" href="#biblio-25"><span class="bib-index">25</span></a>, <a class="bibliocite" name="cite-23" href="#biblio-23"><span class="bib-index">23</span></a>]</span>, also see Eq. (<a class="Reference" href="#himEexact">5↑</a>).
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-5">5</a> quenching the interparticle interaction: MCTDHB and FCI vs. exact results
</h1>
<div class="Unindented">
<a class="Label" name="HIMquench"> </a>
</div>
<div class="Indented">
In the previous section we have seen that the numerically exact ground state solutions of the HIM can be obtained using the MCTDHB method with just a few self-consistent orbitals. The standard full configuration interaction method utilizing the &ldquo;non-optimal&rdquo; fixed-shape orbitals of the non-interacting system has demonstrated a much worser convergence. It makes the usage of the direct diagonalization method with &ldquo;non-optimal&rdquo; orbitals for large particle numbers impractical. Within the number of orbitals technically allowed to be used (this number defines the size of the respective secular matrix to be diagonalized) the quality of the obtained many-body results is unsatisfactory. They can be worser than the one-orbital self-consistent mean-field (GP) results, see Fig. <a class="Reference" href="#HIM_MCTDHB">1↓</a>.
</div>
<div class="Indented">
Having established the great relevance of self-consistency for statics, in the present section we clarify its impact on the quantum many-boson dynamics. The main difference between statics and dynamics is that quantum dynamics involves a lot of excited states and, therefore, the applied many-body method has to be capable to describe them. Indeed, the evolution of any given initial many-body state is obtained as a solution of the time-dependent Schrödinger equation: <div class="formula">
<a class="eqnumber" name="TDEXPAND">(7) </a><span class="environment"><span class="arrayrow">
<span class="arraycell align-r">
Ψ(<i>t</i>)
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="limits"><sup class="limit">∞</sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>j</i> = 0</sub></span><i>a</i><sub><i>j</i></sub><i>e</i><sup> − <i>iE</i><sub><i>j</i></sub><i>t</i></sup>Φ<sub><i>j</i></sub>(<i>X⃗</i>)
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">

</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
<span class="limits"><sup class="limit">∞</sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>j</i> = 0</sub></span>⟨Ψ<sub>0</sub>(<i>X⃗</i>)|Φ<sub><i>j</i></sub>(<i>X⃗</i>)⟩<i>e</i><sup> − <i>iE</i><sub><i>j</i></sub><i>t</i></sup>Φ<sub><i>j</i></sub>(<i>X⃗</i>).
</span>

</span>
</span>
</div>
Here <span class="formula">Ψ<sub>0</sub>(<i>X⃗</i>)</span> is the initial state and <span class="formula">Φ<sub><i>j</i></sub>(<i>X⃗</i>)</span> and <span class="formula"><i>E</i><sub><i>j</i></sub></span> are the eigenstates and respective eigenenergies of the quantum system, <span class="formula"><i>X⃗</i></span> are the coordinates of the constituting particles. For the HIM considered here all the eigenstates and respective eigenenergies are known in the center of mass frame, so, to study the exact evolution of the many-body system one needs to evaluate the overlap integrals <span class="formula"><i>a</i><sub><i>j</i></sub> = ⟨Ψ<sub>0</sub>(<i>X⃗</i>)|Φ<sub><i>j</i></sub>(<i>X⃗</i>)⟩</span>. When computing with the MCTDHB we, of course, work in the laboratory frame and the time-dependent many-body wavefunction is a complicated non-terminating expansion in terms of permanents.
</div>
<div class="Indented">
Let us study a scenario with the HIM Hamiltonian where the many-body dynamics is activated by a sudden quench of the interparticle interaction strength. It is worthwhile to mention that the MCTDHB method has been successfully used in Ref. <span class="bibcites">[<a class="bibliocite" name="cite-31" href="#biblio-31"><span class="bib-index">31</span></a>]</span> to describe such a scenario for ultra-cold systems with contact interaction. On the experimental side the quench of the interparticle interaction is a routine procedure controlled by the Feshbach resonance technique. We assume that the initial state just before the quench was the ground state of the non-interacting system. What kind of dynamics is anticipated in this case? The initial state, i.e., the ground state of the harmonically trapped system is symmetric, implying that the one-body density has &ldquo;gerade&rdquo; symmetry. The sudden quench of the interparticle interaction cannot break this symmetry. So, we expect that a change (quench) of the interparticle interaction leads to a &ldquo;breathing&rdquo; dynamics of the system — the many-body wavefunction changes its shape such that the &ldquo;gerade&rdquo; symmetry is preserved. This dynamical behavior is general and persists in other many-body systems with symmetric trap potentials as well, e.g., in ultra-cold systems with contact interactions.
</div>
<h2 class="Subsection">
<a class="toc" name="toc-Subsection-5.1">5.1</a> Breathing dynamics for <span class="formula"><i>N</i> = 2</span> bosons
</h2>
<div class="Unindented">
We first consider the two-boson HIM system where the interparticle interaction strength <span class="formula"><i>K</i><sub>0</sub></span> is suddenly quenched from zero to <span class="formula"><i>K</i><sub>0</sub> = 0.5</span>. The exact general eigenstates of the HIM system are products of two HO wavefunctions — one describes the motion of the relative <span class="formula"><i>q</i><sub>1</sub> = <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>x</i><sub>2</sub> − <i>x</i><sub>1</sub></span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root">2</span><span class="ignored">)</span></span></span><span class="ignored">)</span></span></span> coordinate, another the motion of the center of mass <span class="formula"><i>q</i><sub>2</sub> = <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>x</i><sub>2</sub> + <i>x</i><sub>1</sub></span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root">2</span><span class="ignored">)</span></span></span><span class="ignored">)</span></span></span>. Clearly, the center of mass part always has bosonic symmetry — it does not change sign when coordinates of the particles are permuted. In contrast, the relative part can be either bosonic (symmetric) or fermionic (antisymmetric), depending on the parity of the Hermite polynomial involved. For example, the first excited HO state of the relative part <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">1</sub></span>(<i>q</i><sub>1</sub>, <i>ω̃</i>) ~ <i>H</i><sub>1</sub>(<span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω̃</i></span><span class="ignored">)</span></span><i>q</i><sub>1</sub>)<i>e</i><sup> − <i>ω̃</i><span class="fraction"><span class="ignored">(</span><span class="numerator"><i>q</i><span class="scripts"><sup class="script">2</sup><sub class="script">1</sub></span></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span></sup> ≡ <span class="scriptfont"> N</span><i>q</i><sub>1</sub><i>e</i><sup> − <i>ω̃</i><span class="fraction"><span class="ignored">(</span><span class="numerator"><i>q</i><span class="scripts"><sup class="script">2</sup><sub class="script">1</sub></span></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span></sup></span> [<span class="formula"><i>ω̃</i> = <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i> + <i>δ</i><sub>2</sub></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span></span>] is fermionic, because the permutation <span class="formula"><i>x</i><sub>1</sub>↔<i>x</i><sub>2</sub></span> changes the sign of the <span class="formula"><i>q</i><sub>1</sub> = <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>x</i><sub>2</sub> − <i>x</i><sub>1</sub></span><span class="ignored">)/(</span><span class="denominator"><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root">2</span><span class="ignored">)</span></span></span><span class="ignored">)</span></span></span> (first Hermite polynomial) and, therefore, the overall sign of the total wavefunction <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">1</sub></span>(<i>q</i><sub>1</sub>, <i>ω̃</i>)<i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script">0</sub></span>(<i>q</i><sub>2</sub>, <i>ω</i>)</span>. Using this argumentation one can conclude that all even excited HO states <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>i</i></sub></span>(<i>q</i><sub>1</sub>, <i>ω̃</i>)</span> with <span class="formula"><i>i</i> = 0, 2, 4, ...</span> of the relative part are bosonic and all odd ones, i.e., <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>i</i></sub></span>(<i>q</i><sub>1</sub>, <i>ω̃</i>)</span> <span class="formula"><i>i</i> = 1, 3, 5, ...</span> are fermionic.
</div>
<div class="Indented">
Having understood the nature of the bosonic and fermionic solutions of the two-particle HIM problem we are ready to analyze the excitations responsible for the &ldquo;breathing&rdquo; dynamics. The &ldquo;ungerade bosonic excitations&rdquo; of the HIM model are activated by odd excitations of the center of mass part <span class="formula"><i>ψ</i><span class="scripts"><sup class="script"><i>ho</i></sup><sub class="script"><i>i</i></sub></span>(<i>q</i><sub>2</sub>, <i>ω</i>)</span>, which oscillates with the original trap frequency <span class="formula"><i>ω</i></span>. The lowest &ldquo;gerade&rdquo; excitation corresponds to the second excited state of the relative part, the next &ldquo;gerade&rdquo; bosonic excitation — to the fourth excited state and so on. In principle, the next class of the &ldquo;gerade&rdquo; excitations appears by product of the HO solutions corresponding to the second excited state of the center of mass motion and every even excitation of the relative motion. In the studied quench dynamics we start from the ground state of the non-interacting system, implying that the center of mass motion is in the ground state. Therefore, all excited states for which the center of mass is excited are orthogonal to such an initial state and, hence, do not contribute to the dynamics; this is because the overlap integrals of these states with the initial state are zero. Summarizing, the sudden quench of the interparticle interaction in the HIM leads to the breathing dynamics with breathing frequencies <span class="formula"><i>ω</i>(<i>n</i>) = 2<i>n</i><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><sup>2</sup> + 4<i>K</i><sub>0</sub></span><span class="ignored">)</span></span></span>, with main excitation frequency <span class="formula"><i>ω</i>(<i>n</i> = 1) ≡ <i>ω</i><sub><i>breath</i></sub></span> and all its overtones with <span class="formula"><i>n</i> = 2, 3, 4, ...</span>. These frequencies are obtained as energy differences between the ground and respective excited eigenstates.
</div>
<div class="Indented">
We use the Mathematica package <span class="bibcites">[<a class="bibliocite" name="cite-43" href="#biblio-43"><span class="bib-index">43</span></a>]</span> to compute the required overlap integrals in Eq. (<a class="Reference" href="#TDEXPAND">7↑</a>) and to get the exact time-dependent two-boson wavefunction. Here we have to mention that instead of the infinite summation the contributions from the 60 exact lowest-in-energy excited states are taken into account; this is more than sufficient for numerical convergence. Next, the exact one-body density as a function of time is obtained according to its definition: the two-body wavefunction is multiplied by its complex conjugate and one coordinate is integrated out. In Fig. <a class="Reference" href="#HIM_QUENCH_N2">2↓</a> we plot the exact value of the density at the trap center as a function of time by a bold red line. The first two breathing cycles are depicted in the left panel of this figure, the right panel presents the breathing dynamics at longer propagation times. The numerical <span class="formula"><i>M</i></span>-orbital MCTDHB and FCI results are depicted by bold symbols and dotted lines, respectively.
</div>
<div class="Indented">
The first observation seen from Fig. <a class="Reference" href="#HIM_QUENCH_N2">2↓</a> is that the exact density at the middle of the trap oscillates periodically with the breathing frequency <span class="formula"><i>ω</i><sub><i>breath</i></sub> = 2<span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><sup>2</sup> + 4<i>K</i><sub>0</sub></span><span class="ignored">)</span></span></span>. However, the shape of the oscillation differs from the simple <span class="formula"> ~ cos(<i>ω</i><sub><i>breath</i></sub><i>t</i>)</span> function plotted to guide the eye by a solid black line. This is the result of the contributions from the overtones originating from the higher excited states. The two-orbital MCTDHB(2) solution provides essentially an exact description of the dynamics till half of the breathing cycle — notice the triangles following the exact results. The three-orbital MCTDHB(3) results, plotted by filled circles, are on-top of the exact curve for the first breathing cycle; small deviations from the exact results become visible at the second breathing cycle. The MCTDHB(4) with four time-adaptive orbitals gives the exact description of the first two breathing cycles. The FCI dynamics with four fixed-shape orbitals, plotted by a double-dashed line, starts to deviate from the exact results already at very short times. Even the six-orbital FCI dynamics, depicted by a dashed line, starts to deviate from the exact result after one third of the first breathing cycle. A quite accurate description of the first two breathing oscillations is only obtained on the eight-orbital FCI level. Summarizing, to describe the first two breathing oscillations of the two-particle HIM problem one needs either four time-adaptive orbitals [MCTDHB(4)] or eight fixed-shape orbitals [FCI(8)]. The above analysis also shows that with time more time-adaptive orbitals are needed to describe the exact dynamics.
</div>
<div class="Indented">
Now we quantify the performance of the MCTDHB and FCI methods to describe the quantum dynamics at longer times. In the right part of Fig. <a class="Reference" href="#HIM_QUENCH_N2">2↓</a> we plot the oscillations of the density at the middle of the trap at longer times. Exact results are depicted by a bold line, the MCTDHB results are depicted by filled symbols and the FCI results are plotted by dashed lines. The exact density continues to oscillate with the main <span class="formula"><i>ω</i><sub><i>breath</i></sub></span> and its overtones. The MCTDHB(<span class="formula"><i>M</i> &gt; 5</span>) results are numerically exact. The FCI(12) result, plotted by a dense-dashed line, follows the exact one, while the FCI(6) and FCI(7) are clearly off.
</div>
<div class="Indented">
Concluding, for longer propagation time one has to span a larger Fock subspace, i.e., one has to use a larger number of one-particle functions to build the permanents. The difference between the quantum dynamics utilizing fixed-shape and time-adaptive orbital basis sets used is clearly seen — to gain a desired accuracy of the propagation one needs to use at least twice as many fixed-shape orbitals than time-adaptive ones. The second important observation is that if a desired convergence of the many-body dynamics is achieved on the <span class="formula"><i>M</i></span>-orbital level, further extension of the Fock space is unnecessary, the inclusion of the extra orbitals does not impact the result. This is a general consequence of the variational principle used. It is known for the FCI method and now proven for the MCTDHB method, which is based on the time-dependent Dirac-Frenkel variational principle. This feature allows us to define a practical strategy for MCTDHB computations. If the dynamics obtained on the MCTDHB(<span class="formula"><i>M</i></span>) and MCTDHB(<span class="formula"><i>M</i> + 1</span>) levels are identical we conclude that numerical convergence to the exact results is reached. In other words, the many-body wavefunction built from <span class="formula"><i>M</i></span> time-adaptive orbitals is the converged solution of the time-dependent many-boson Schrödinger equation.
</div>
<h2 class="Subsection">
<a class="toc" name="toc-Subsection-5.2">5.2</a> Breathing dynamics for <span class="formula"><i>N</i> = 10</span> bosons
</h2>
<div class="Unindented">
Now we examine and compare the performance of the MCTDHB and FCI methods to treat the time-dependent dynamics of the HIM system with N=10 bosons for the same quench scenario as studied before for a system with N=2 bosons. By analyzing the structure of the excited states for the <span class="formula"><i>N</i> = 10</span> system we arrive at the conclusion that a sudden quench leads to the many-body breathing dynamics with frequencies <span class="formula"><i>ω</i>(<i>n</i>) = 2<i>n</i><span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><sup>2</sup> + 2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span></span>, obtained as energy differences between the ground and respective excited eigenstates. Here, the lowest excitation <span class="formula"><i>ω</i>(<i>n</i> = 1) ≡ <i>ω</i><sub><i>breath</i></sub></span> is responsible for the main breathing excitation frequency, higher excited states with <span class="formula"><i>n</i> = 2, 3, 4, ...</span> result in overtones. The exact results are in principle available for the systems with any number of particles <span class="bibcites">[<a class="bibliocite" name="cite-25" href="#biblio-25"><span class="bib-index">25</span></a>, <a class="bibliocite" name="cite-23" href="#biblio-23"><span class="bib-index">23</span></a>]</span>. However the straightforward way of evaluation of the exact time-dependent many-body wavefunction and the respective density successfully applied for two-boson is much more involved for the ten-particle system. Hence, we employ the numerical MCTDHB method. The first three breathing cycles are depicted in the left part of Fig. <a class="Reference" href="#HIM_QUENCH_N10">3↓</a>. The MCTDHB results utilizing <span class="formula"><i>M</i> = 8, 9, 10</span> time-adaptive orbitals plotted by dashed lines are indistinguishable from each other. The computational strategy verified and proven above allows us to conclude that the numerically exact description of the TDSE is achieved.
</div>
<div class="Indented">
Another important observation seen in Fig. <a class="Reference" href="#HIM_QUENCH_N10">3↓</a> is that the numerically exact MCTDHB results deviate substantially from a simply fitted <span class="formula"> ~ cos(<i>ω</i><sub><i>breath</i></sub><i>t</i>)</span> curve, plotted by a solid black line to guide the eye. This is direct evidence that the contribution from higher excited states, responsible for higher overtones, to the breathing dynamics of the <span class="formula"><i>N</i> = 10</span> boson system is much stronger than it was in the <span class="formula"><i>N</i> = 2</span> system studied before [compare the exact curve and fit to the <span class="formula"> ~ cos(<i>ω</i><sub><i>breath</i></sub><i>t</i>)</span> curve in Fig. <a class="Reference" href="#HIM_QUENCH_N2">2↓</a>]. The FCI result with <span class="formula"><i>M</i> = 16</span> fixed-shape orbitals depicted by filled triangles follows the numerically exact MCTDHB curves only for a very short initial time — for about one half of the first breathing cycle. For longer propagation times the FCI(16) predictions deviate from the exact results.
</div>
<div class="Indented">
The right part of Fig. <a class="Reference" href="#HIM_QUENCH_N10">3↓</a> depicts on an enlarged scale the breathing dynamics at longer times. At all plotted propagation times the eight-orbitals MCTDHB(8) method provides a very accurate description of the many-boson dynamics while the MCTDHB computations with <span class="formula"><i>M</i> = 9, 10</span> time-adaptive orbitals are numerically exact. We conclude that the usage of time-adaptive orbitals provides an enormous benefit for the accurate description of quantum dynamics of systems with larger particle numbers. In contrast to the MCTDHB method which is capable to provide numerically exact results with a few time-adaptive orbitals, the FCI treatments even with much larger Fock subspaces spanned can not provide even a qualitative description of the dynamics.
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-6">6</a> The time-dependent HIM: Non-equilibrium dynamics
</h1>
<div class="Unindented">
<a class="Label" name="TDHIM"> </a>
</div>
<div class="Indented">
The above discussed visualization of the HIM system as a medium made of <span class="formula"><i>N</i> − 1</span> non-interacting &ldquo;relative&rdquo; particles in which the effective mass particle (representing the center of mass coordinate) lives, allows for a simple physical time-dependent generalization. Without loss of the separability one can assume that the effective-mass particle moves now not in the stationary but in a time-dependent harmonic potential with driving frequency <span class="formula"><i>ω</i>(<i>t</i>)</span>. Moreover, we assume that during this motion the medium representing the relative coordinates remains undisturbed <span class="formula"><i>δ</i><sub><i>N</i></sub> = <i>const</i>.</span>. Surprisingly, the Hamiltonian corresponding to this problem takes on a simple form in the center of mass frame: <div class="formula">
<a class="eqnumber" name="omtranstd1">(8) </a><i>Ĥ</i>(<i>t</i>) = <i>Ĥ</i><sub><i>rel</i></sub> + <i>Ĥ</i><sub><i>CM</i></sub>(<i>t</i>) = <span class="limits"><sup class="limit"><i>N</i> − 1</sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> = 1</sub></span>( − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>q⃗</i><sub><i>i</i></sub></sub></span> + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>δ</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>N</i></sub></span><i>q⃗</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>i</i></sub></span>) − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>q⃗</i><sub><i>N</i></sub></sub></span> + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i><sup>2</sup>(<i>t</i>)<i>q⃗</i><span class="scripts"><sup class="script">2</sup><sub class="script"><i>N</i></sub></span>.
</div>
One can apply a reverse engineering and transform this new time-dependent HIM problem back to the laboratory frame: <div class="formula">
<a class="eqnumber" name="omtranstd2">(9) </a><i>Ĥ</i>(<i>t</i>) = <span class="limits"><sup class="limit"><i>N</i></sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> = 1</sub></span>[ − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>r<sub>i</sub>⃗</i></sub></span> + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i>(<i>t</i>)<sup>2</sup><i>r<sub>i</sub>⃗</i><sup>2</sup>] + <i>K</i>(<i>t</i>)<span class="limits"><sup class="limit"><i>N</i></sup><span class="limit">⎲</span><span class="limit">⎳</span><sub class="limit"><i>i</i> &lt; <i>j</i></sub></span><span class="symbol">(</span><i>r⃗</i><sub><i>i</i></sub> − <i>r⃗</i><sub><i>j</i></sub><span class="symbol">)</span><sup>2</sup>.
</div>
This coupled time-dependent Hamiltonian corresponds to the situation where all &ldquo;real&rdquo; particles are trapped in the time-dependent potential <span class="formula"><i>V̂</i>(<i>r⃗</i>, <i>t</i>) = <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i>(<i>t</i>)<sup>2</sup><i>r⃗</i><sup>2</sup></span> and interact via time-dependent harmonic interparticle interaction potential of strength <span class="formula"><i>K</i>(<i>t</i>)</span> [which depends on <span class="formula"><i>ω</i>(<i>t</i>)</span> and <span class="formula"><i>δ</i><sub><i>N</i></sub></span>]. For the external trapping potential driven by a time-dependent function <span class="formula"><i>f</i>(<i>t</i>)</span>: <div class="formula">
<a class="eqnumber" name="omtrans">(10) </a><i>ω</i>(<i>t</i>) = <i>ω</i><sub>0</sub><span class="symbol">[</span>1 + <i>f</i>(<i>t</i>)<span class="symbol">]</span>.
</div>
the imposed above requirement <span class="formula"><i>δ</i><sub><i>N</i></sub> = <span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><span class="scripts"><sup class="script">2</sup><sub class="script">0</sub></span> + 2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span> = <i>const</i>.</span> implies that the interparticle interaction strength has to be driven with the &ldquo;compensating&rdquo; time-dependency: <div class="formula">
<a class="eqnumber" name="Ktrans">(11) </a><i>K</i>(<i>t</i>) = <i>K</i><sub>0</sub><span class="array"><span class="arrayrow"><span class="bracket align-left">⎡</span></span><span class="arrayrow"><span class="bracket align-left">⎣</span></span></span>1 − <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i><span class="scripts"><sup class="script">2</sup><sub class="script">0</sub></span></span><span class="ignored">)/(</span><span class="denominator">2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span><i>f</i>(<i>t</i>)<span class="array"><span class="arrayrow"><span class="bracket align-right">⎤</span></span><span class="arrayrow"><span class="bracket align-right">⎦</span></span></span>.
</div>
Since the Hamiltonian (<a class="Reference" href="#omtranstd1">8↑</a>) or (<a class="Reference" href="#omtranstd2">9↑</a>) is now time-dependent, the total energy is, of course, no longer conserved.
</div>
<div class="Indented">
Let us consider a situation where the medium representing <span class="formula"><i>N</i> − 1</span> relative particles is in the ground state of the harmonic potential with frequency <span class="formula"><i>δ</i><sub><i>N</i></sub> = <span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><span class="scripts"><sup class="script">2</sup><sub class="script">0</sub></span> + 2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span></span>. Its energy is the time-independent constant <span class="formula"><span class="fraction"><span class="ignored">(</span><span class="numerator"><i>D</i></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>N</i> − 1)<i>δ</i><sub><i>N</i></sub></span>. The time-dependency of the full problem then originates from the driving of the center of mass: <div class="formula">
<a class="eqnumber" name="TDSE1p">(12) </a><i>Ĥ</i><sub><i>CM</i></sub><i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>) =  − <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>∂<span class="scripts"><sup class="script">2</sup><sub class="script"><i>q<sub>N</sub>⃗</i></sub></span><i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>) + <span class="fraction"><span class="ignored">(</span><span class="numerator">1</span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span><i>ω</i><span class="scripts"><sup class="script">2</sup><sub class="script">0</sub></span>[1 + <i>f</i>(<i>t</i>)]<sup>2</sup><i>q<sub>N</sub>⃗</i><sup>2</sup><i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>) = <i>i</i><span class="fraction"><span class="ignored">(</span><span class="numerator">∂</span><span class="ignored">)/(</span><span class="denominator">∂<i>t</i></span><span class="ignored">)</span></span><i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>).
</div>
The solution <span class="formula"><i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>)</span> of this <i>one-particle</i> Schrödinger equation can easily be obtained numerically, see Refs. <span class="bibcites">[<a class="bibliocite" name="cite-24" href="#biblio-24"><span class="bib-index">24</span></a>, <a class="bibliocite" name="cite-41" href="#biblio-41"><span class="bib-index">41</span></a>]</span>. The final expression for the expectation value of the Hamiltonian Eq. (<a class="Reference" href="#omtranstd1">8↑</a>) or  (<a class="Reference" href="#omtranstd2">9↑</a>) reads: <div class="formula">
<a class="eqnumber" name="TDEXPVALE">(13) </a>⟨Ψ(<i>t</i>)∣<i>Ĥ</i><sub><i>rel</i></sub> + <i>Ĥ</i><sub><i>CM</i></sub>(<i>t</i>)∣Ψ(<i>t</i>)⟩ = <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>D</i></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>N</i> − 1)<i>δ</i><sub><i>N</i></sub> + <i>ϵ</i>(<i>t</i>), 
</div>
where <span class="formula"><i>ϵ</i>(<i>t</i>) = ⟨<i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>)|<i>Ĥ</i><sub><i>CM</i></sub>(<i>t</i>)|<i>ψ</i>(<i>q<sub>N</sub>⃗</i>, <i>t</i>)⟩</span>. Interestingly, the special kind of time-dependency used in Eqs. (<a class="Reference" href="#omtrans">10↑</a>, <a class="Reference" href="#Ktrans">11↑</a>) also implies that the time-dependent part <span class="formula"><i>ϵ</i>(<i>t</i>)</span> of the expectation value of the total Hamiltonian <span class="formula"><i>Ĥ</i>(<i>t</i>)</span> depends neither on the number of particles <span class="formula"><i>N</i></span> nor on the interaction strength <span class="formula"><i>K</i><sub>0</sub></span>. So, the systems with different particle numbers <span class="formula"><i>N</i></span> and different interparticle interaction strengths <span class="formula"><i>K</i><sub>0</sub></span> possess the same time-dependent fraction.
</div>
<div class="Indented">
It is instructive to state here that the time-dependencies, Eqs. (<a class="Reference" href="#omtrans">10↑</a>,<a class="Reference" href="#Ktrans">11↑</a>), can be more general and it is of course not necessary to choose them such that the relative Hamiltonian <span class="formula"><i>Ĥ</i><sub><i>rel</i></sub></span> remains time-independent, i.e., keeping <span class="formula"><i>δ</i><sub><i>N</i></sub> = <i>const</i>.</span> Yet, there is one important advantage to this particular choice: in the center of mass frame, Eq. (<a class="Reference" href="#omtranstd1">8↑</a>), the relative part is known analytically and to solve the problem completely and exactly one needs to integrate only a single one-particle Schrödinger equation, Eq. (<a class="Reference" href="#TDSE1p">12↑</a>). In contrast, to find the solution in the laboratory frame one has to solve the time-dependent many-boson Schrödinger equation with a <i>time-dependent trap potential</i> and <i>time-dependent interparticle interactions</i>. While the former task is a manageable standard routine, the latter one comprises a very involved and appealing theoretical and numerical problem. The main goal of the present section is to show that the MCTDHB method is capable to tackle time-dependent scenarios numerically exactly even in the most involved setups: time-dependent one-particle potentials <span class="formula"><i>V̂</i>(<i>r⃗</i>, <i>t</i>)</span> and time-dependent two-body interactions <span class="formula"><i>Ŵ</i>(<i>r⃗</i>, <i>r⃗</i>’, <i>t</i>)</span>.
</div>
<div class="Indented">
In what follows, we investigate the dynamics of the one-dimensional HIM system with time-dependent trapping (<a class="Reference" href="#omtrans">10↑</a>) and interaction (<a class="Reference" href="#Ktrans">11↑</a>) potentials driven by two different functions <div class="formula">
<a class="eqnumber" name="f1_f2">(14) </a><span class="environment"><span class="arrayrow">
<span class="arraycell align-r">
<i>f</i><sub>1</sub>(<i>t</i>)
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
0.2sin<sup>2</sup>(<i>t</i>), 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
 
</span>
<span class="arraycell align-c">
 
</span>
<span class="arraycell align-l">
 
</span>

</span>
<span class="arrayrow">
<span class="arraycell align-r">
<i>f</i><sub>2</sub>(<i>t</i>)
</span>
<span class="arraycell align-c">
 = 
</span>
<span class="arraycell align-l">
sin(<i>t</i>)cos(2<i>t</i>)sin(0.5<i>t</i>)sin(0.4<i>t</i>).
</span>

</span>
</span>
</div>
The one-body center-of-mass Schrödinger equation (<a class="Reference" href="#TDSE1p">12↑</a>) with the respective time-dependent potentials is integrated numerically exactly to obtain the corresponding one-body energies <span class="formula"><i>ϵ</i><sub>1</sub>(<i>t</i>), <i>ϵ</i><sub>2</sub>(<i>t</i>)</span>.
</div>
<div class="Indented">
Let us first study the time-dependent HIM system made of <span class="formula"><i>N</i> = 10</span> bosons with a relatively simple periodic driving function <span class="formula"><i>f</i><sub>1</sub>(<i>t</i>)</span> and <span class="formula"><i>K</i><sub>0</sub> = 0.5</span>. In the lower part of Fig. <a class="Reference" href="#HIMTD_ft1">4↓</a> we plot the time-dependent part <span class="formula"><i>ϵ</i><sub>1</sub>(<i>t</i>)</span> of the respective expectation value of the total Hamiltonian <span class="formula"><i>Ĥ</i>(<i>t</i>)</span> computed by using different levels of the MCTDHB(<span class="formula"><i>M</i></span>) theory. The numerically exact results for <span class="formula"><i>ϵ</i><sub>1</sub>(<i>t</i>)</span> depicted by open circles are obtained by solving directly the one-particle time-dependent Schrödinger equation. It is important to notice that the oscillatory motion of the center of mass results in a relatively small contribution to the total energy: the value of <span class="formula"><i>ϵ</i><sub>1</sub>(<i>t</i>)</span> is of the order of a single-particle energy.
</div>
<div class="Indented">
The close inspection of Fig. <a class="Reference" href="#HIMTD_ft1">4↓</a> shows that the MCTDHB(3) computation, depicted by a dashed-doted line, follows the exact curve until <span class="formula"><i>t</i> ≈ 5</span>. To describe the exact dynamics for longer propagation times one needs to use higher levels of the MCTDHB(<span class="formula"><i>M</i></span>) theory, i.e., more time-adaptive orbitals are needed. The MCTDHB(5) result, plotted by a dense-dashed line, is exact until <span class="formula"><i>t</i> ≈ 15</span>, while the MCTDHB(6) one depicted by a simple dashed line is exact until <span class="formula"><i>t</i> ≈ 30</span>. The double-dashed line depicting the MCTDHB(7) result reproduces the exact time-dependency of the total energy at all the times considered here.
</div>
<div class="Indented">
In the context of ultra-cold physics the Gross-Pitaevskii mean-field theory is considered as one of the main working tool to describe the dynamics of bosonic systems with time-dependent traps and time-dependent interactions. In Fig. <a class="Reference" href="#HIMTD_ft1">4↓</a> we also plot the results obtained by solving GP equation with the finite-range harmonic interaction, which is identical to the lowest level MCTDHB(1) theory, by a bold solid line. The GP theory is incapable to describe the time-dependent energy correctly even for short times. Note that <span class="formula"><i>N</i> = 10</span> only.
</div>
<div class="Indented">
Now we study the HIM systems made of N=10 and N=50 bosons with <span class="formula"><i>K</i><sub>0</sub> = 0.5</span> driven by quite a complicated time-dependent function <span class="formula"><i>f</i><sub>2</sub>(<i>t</i>)</span>, depicted in the upper part of Fig. <a class="Reference" href="#HIMTD_ft2">5↓</a>. The exact <span class="formula"><i>ϵ</i><sub>2</sub>(<i>t</i>)</span>, obtained by solving the corresponding one-particle Schrödinger equation, is plotted by open red circles. As it was discussed above, this time-dependent fraction of the expectation value of the Hamiltonian <span class="formula"><i>Ĥ</i>(<i>t</i>)</span> is the same for both systems. We use different levels of the MCTDHB(<span class="formula"><i>M</i></span>) theory to compute the time-dependent contribution <span class="formula"><i>ϵ</i><sub>2</sub>(<i>t</i>)</span> to the energy of the studied systems. The dashed-dotted and dashed-double-dotted lines are used, respectively, to depict the MCTDHB(<span class="formula"><i>M</i> = 6</span>) and MCTDHB(<span class="formula"><i>M</i> = 7</span>) results for the system with <span class="formula"><i>N</i> = 10</span> bosons. The MCTDHB(<span class="formula"><i>M</i> = 5</span>) and MCTDHB(<span class="formula"><i>M</i> = 6</span>) results for the system made of <span class="formula"><i>N</i> = 50</span> bosons are depicted by the dashed and dense-dashed lines, correspondingly. All the presented numerical MCTDHB results with <span class="formula"><i>M</i> &gt; 5</span> follow the exact lines till <span class="formula"><i>t</i> ≈ 25</span>, indicating that numerical convergence is reached. In Fig. <a class="Reference" href="#HIMTD_ft2">5↓</a> we also depict the corresponding Gross-Pitaevskii results. The GP theory, usually considered to be applicable for systems made of a larger number of particles, provides a semi-qualitative description of the very short initial dynamics (<span class="formula"><i>t</i> ≈ 1</span>), afterwards its predictions sharply deteriorates.
</div>
<div class="Indented">
Summarizing, the MCTDHB(<span class="formula"><i>M</i></span>) computations with a <i>given</i> number of time-adaptive orbitals start to deviate from the exact result with time, see Fig. <a class="Reference" href="#HIMTD_ft2">5↓</a> and its inset. The time-dependent variational principle used in the MCTDHB method implies that the MCTDHB computations done with a larger number of the time-adaptive orbitals remain &ldquo;on-top&rdquo; of the exact curve for <i>longer</i> propagation times. Even when the numerical many-body results slightly deviate form the exact at longer propagation times they are quantitative and quite accurate -- all the spectral features of the exact behavior are reproduced, see Fig. <a class="Reference" href="#HIMTD_ft2">5↓</a>. In conclusion, the MCTDHB method is capable to provide numerically converged results for time-dependent Hamiltonians with very general driving scenarios, where both the external trap and interparticle interactions are driven in quite a complicated way.
</div>
<h1 class="Section">
<a class="toc" name="toc-Section-7">7</a> Discussion and Outlook
</h1>
<div class="Unindented">
<a class="Label" name="sum"> </a>
</div>
<div class="Indented">
In the present work we have compared the quantum many-boson physics of the HIM system described in the laboratory and in the center of mass frames. In contrast to the center of mass frame where the HIM problem is exactly solvable, in the laboratory frame one has to apply numerical efforts to reproduce the exact results, even in the two particle case. The relevance of self-consistency and time-adaptivity is demonstrated. To solve time-independent problems, the standard many-body full configuration interaction (exact diagonalization) method requires to use a large number of fixed-shape (non-optimal) one-particle basis functions, thereby restricting its applicability to few-particle systems. The usage of the MCTDHB method utilizing variational self-consistent basis set allows to attack systems with larger particle numbers. To verify how good the obtained static solution is one can use the straightforward methodology: by comparing the eigenstates obtained by imaginary time propagation of the MCTDHB equations with M and M+1 orbitals one can conclude that numerical convergence is achieved. See table <a class="Reference" href="#tabHIMenergies">1↑</a> for reference.
</div>
<div class="Indented">
To check the relevance of the time-adaptivity we have first studied the time-dependent HIM problem where the dynamics are initiated by a sudden quench of the interparticle interaction strength from zero to some finite value. It has been shown that the many-body method (MCTDHB) utilizing variationally optimal time-adaptive orbitals allows to obtain the numerically exact solutions for much longer propagation times in comparison to the fixed-orbital FCI method spanning a Fock space of the same size. The methodology in determining the accuracy of the time-dependent solution obtained is to compare the properties of the MCTDHB solutions computed by using M and M+1 time-adaptive orbitals at different propagation times. If more time-adaptive orbitals are used than needed, the Dirac-Frenkel variational principle keeps the superfluous orbitals unoccupied, i.e., they do not contribute to the now converged and exact many-boson wavefunction. Generally, we have found that one needs less time-adaptive orbitals to converge the results for bosonic systems with increasing number <span class="formula"><i>N</i></span> of particles when the interaction strength <span class="formula">Λ = <i>K</i><sub>0</sub>(<i>N</i> − 1)</span> is kept fixed.
</div>
<div class="Indented">
In the broader context several other methods of the family of the multiconfigurational methods have been benchmarked in the field. The first of these, the multiconfigurational time-dependent Hartree (MCTDH, MCTDHB’s mother method) <span class="bibcites">[<a class="bibliocite" name="cite-28" href="#biblio-28"><span class="bib-index">28</span></a>, <a class="bibliocite" name="cite-39" href="#biblio-39"><span class="bib-index">39</span></a>, <a class="bibliocite" name="cite-1" href="#biblio-1"><span class="bib-index">1</span></a>]</span>, was benchmarked with standard wave-packet propagation <span class="bibcites">[<a class="bibliocite" name="cite-28" href="#biblio-28"><span class="bib-index">28</span></a>, <a class="bibliocite" name="cite-39" href="#biblio-39"><span class="bib-index">39</span></a>, <a class="bibliocite" name="cite-1" href="#biblio-1"><span class="bib-index">1</span></a>]</span> as well as with experimental spectra (see, e.g., Refs. <span class="bibcites">[<a class="bibliocite" name="cite-7" href="#biblio-7"><span class="bib-index">7</span></a>, <a class="bibliocite" name="cite-40" href="#biblio-40"><span class="bib-index">40</span></a>]</span>). The MCTDH for fermions (MCTDHF, a sister method of MCTDHB) <span class="bibcites">[<a class="bibliocite" name="cite-22" href="#biblio-22"><span class="bib-index">22</span></a>, <a class="bibliocite" name="cite-20" href="#biblio-20"><span class="bib-index">20</span></a>, <a class="bibliocite" name="cite-38" href="#biblio-38"><span class="bib-index">38</span></a>, <a class="bibliocite" name="cite-27" href="#biblio-27"><span class="bib-index">27</span></a>]</span> was benchmarked with direct numerical solutions of the Schrödinger equation (see, e.g., Refs. <span class="bibcites">[<a class="bibliocite" name="cite-18" href="#biblio-18"><span class="bib-index">18</span></a>, <a class="bibliocite" name="cite-19" href="#biblio-19"><span class="bib-index">19</span></a>]</span>). Similarly, the aim of the present study was to benchmark and assess the properties of the convergence of MCTDHB with respect to the number of variational parameters used. Throughout this work, the MCTDHB method has been benchmarked with the standard HIM. The convergence of the ground state and non-equilibrium dynamics has been demonstrated. We prove, thereby, that the MCTDHB can be used to obtain <i>numerically exact</i> solutions of the many-boson TDSE.
</div>
<div class="Indented">
We have also shown that the exactly-solvable many-body HIM problem can be extended to a driven time-dependent Hamiltonian. Namely, if the time-dependent modulation of the harmonic trap is accompanied by the modulation of the interparticle interaction with the same driving function (with a different amplitude) all the internal excitations in such a system can be compensated and the many-body system behaves as a driven single particle system. For systems with large particle numbers this driving contribution is of the order of a single-particle energy. Physically, it means that the modulation of the harmonic trap can be almost completely compensated by the corresponding modulation of the interparticle interaction.
</div>
<div class="Indented">
The driving scenario proposed for the HIM model is based on the separability of the relative and center of mass coordinates and, therefore, it can also be adapted to other many-body system with such a separability. In particular, it can work in many-body systems trapped in harmonic potentials and interacting via other two-body potentials which depend on the interparticle separation. However it should be noted that, whereas the &ldquo;compensating&rdquo; relation between the trap and interparticle modulations are of simple form for harmonic interactions, see Eqs. (<a class="Reference" href="#omtrans">10↑</a>,<a class="Reference" href="#Ktrans">11↑</a>), in the case of ultra-cold gases (with contact interactions) this relation is expected to be much more involved. Summarizing, in trapped many-particle systems where the center of mass is separable a novel phenomena of &ldquo;dynamical compensation&rdquo; can take place — all the excitations originating from a driven trapping potential can be almost completely dynamically compensated by the respective driving of the interparticle interaction potential.
</div>
<div class="Indented">
The driven many-body system can in principle be realized in the context of ultra-cold physics. It would correspond to an experimental setup where the trap potential (magneto-optical trap) and external magnetic field, used by the Feshbach resonance technique, are driven such that the relative phase and amplitude of the time-dependent modulations can be tuned. By measuring, e.g., the density response as a function of the amplitude of the modulation applied one can scan for and verify the predicted effect. If the dynamical compensation does not take place, more and more excited states would contribute to the dynamics, and in time the density will oscillate with larger and larger amplitude. On the contrary, when the compensation is achieved the density response to the applied modulations remains very weak even for long exposition times. It is important to note that this prediction is valid for many-boson systems where the relative motion is not only in the ground but also in excited states. In particular, it means that the dynamical compensation can work at non-zero temperatures as well.
</div>
<div class="Indented">
<span class="unknown">\acknowledgments</span>
</div>
<div class="Indented">
Computation time on the bwGRiD and the Cray XE6 cluster &ldquo;Hermit&rdquo; at the HLRS in Stuttgart, and financial support by the HGS MathComp, Minerva Short Term research grant, and the DFG also within the framework of the &ldquo;Enable fund&rdquo; of the excellence initiative at Heidelberg university are greatly acknowledged.
</div>
<div class="Indented">
<p><br/>
</p>
<h1 class="biblio">
References
</h1>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-1"><span class="bib-index">1</span></a>] </span> <i><span class="bib-title">Multidimensional Quantum Dynamics: MCTDH Theory and Applications</span></i> (<span class="bib-editor">H.-D. Meyer and F. Gatti and G. A. Worth</span>, ed.). <span class="bib-publisher">WILEY-VCH, Weinheim</span>, <span class="bib-year">2009</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-2"><span class="bib-index">2</span></a>] </span> <span class="bib-authors">J. M. Vogels A. Görlitz, W. Ketterle</span>: “<span class="bib-title">Realization of Bose-Einstein Condensates in Lower Dimensions</span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">130402</span>, <span class="bib-year">2001</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-3"><span class="bib-index">3</span></a>] </span> <span class="bib-authors">J. Werner A. Griesmaier, T. Pfau</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">PRL</span></i>, pp. <span class="bib-pages">160401</span>, <span class="bib-year">2005</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-4"><span class="bib-index">4</span></a>] </span> <span class="bib-authors">A. I. Streltsov, K. Sakmann, A. U. J. Lode, O. E. Alon, L. S. Cederbaum</span>: <i><span class="bib-title">The Multiconfigurational time-dependent Hartree for Bosons package, version 2.1, Heidelberg; http://mctdhb.org</span></i>. <span class="bib-year">2011</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-5"><span class="bib-index">5</span></a>] </span> <span class="bib-authors">O. E. Alon A. I. Streltsov, L. S. Cederbaum</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">130401</span>, <span class="bib-year">2008</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-6"><span class="bib-index">6</span></a>] </span> <span class="bib-authors">A. J. Leggett</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">Rev. Mod. Phys.</span></i>, pp. <span class="bib-pages">307</span>, <span class="bib-year">2001</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-7"><span class="bib-index">7</span></a>] </span> <span class="bib-authors">A. Raab, G. A. Worth, H.-D. Meyer, L. S. Cederbaum</span>: “<span class="bib-title">Molecular dynamics of pyrazine after excitation to the S_2 electronic state using a realistic 24-mode model Hamiltonian</span>”, <i><span class="bib-journal">J. Chem. Phys</span></i>, pp. <span class="bib-pages">936</span>, <span class="bib-year">1999</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-8"><span class="bib-index">8</span></a>] </span> <span class="bib-authors">Iva B ř řezinová, Axel U. J. Lode, Alexej I. Streltsov, Ofir E. Alon, Lorenz S. Cederbaum, Joachim Burgdörfer</span>: “<span class="bib-title">Wave chaos as signature for depletion of a Bose-Einstein condensate</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">013630</span>, <span class="bib-year">2012</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-9"><span class="bib-index">9</span></a>] </span> <span class="bib-authors">C. C. Bradley, C. A. Sackett, J. J. Tollett, R. G. Hulet</span>: “<span class="bib-title">Evidence fof Bose-Einstein Condensation in an Atomic Gas with Attractive Interactions</span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">1687</span>, <span class="bib-year">1995</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-10"><span class="bib-index">10</span></a>] </span> <span class="bib-authors">C. J. Pethick, H. Smith</span>: <i><span class="bib-title">Bose-Einstein condensation in dilute gases</span></i>. <span class="bib-publisher">Cambridge University Press</span>, <span class="bib-year">2008</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-11"><span class="bib-index">11</span></a>] </span> <span class="bib-authors">R. Grimm Ch. Chin, E. Tiesinga</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">RevModPhys</span></i>, pp. <span class="bib-pages">1225</span>, <span class="bib-year">2010</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-12"><span class="bib-index">12</span></a>] </span> <span class="bib-authors">E.A. Cornell, C.E. Wieman</span>: “<span class="bib-title">Nobel Lecture: Bose-Einstein condensation in a dilute gas, the first 70 years and some recent experiments</span>”, <i><span class="bib-journal">Rev. Mod. Phys.</span></i>, pp. <span class="bib-pages">875-893</span>, <span class="bib-year">2002</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-13"><span class="bib-index">13</span></a>] </span> <span class="bib-authors">S. Giorgini F. Dalfovo, S. Stringari</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">Rev. Mod. Phys.</span></i>, pp. <span class="bib-pages">463</span>, <span class="bib-year">1999</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-14"><span class="bib-index">14</span></a>] </span> <span class="bib-authors">M. Greiner, I. Bloch, O. Mandel, T. W. Hänsch, T Esslinger</span>: “<span class="bib-title">Exploring Phase Coherence in a 2D Lattice of Bose-Einstein Condensates</span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">160405</span>, <span class="bib-year">2001</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-15"><span class="bib-index">15</span></a>] </span> <span class="bib-authors">J. Grond, U. Hohenester, J. Schmiedmayer, A. Smerzi</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">023619</span>, <span class="bib-year">2011</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-16"><span class="bib-index">16</span></a>] </span> <span class="bib-authors">Julian Grond, Jörg Schmiedmayer, Ulrich Hohenester</span>: “<span class="bib-title">Optimizing number squeezing when splitting a mesoscopic condensate</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">021603(R)</span>, <span class="bib-year">2009</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-17"><span class="bib-index">17</span></a>] </span> <span class="bib-authors">K. Henderson, C. Ryu, C. MacCormick, M. G. Boshier</span>: “<span class="bib-title">Experimental demonstration of painting arbitrary and dynamic potentials for Bose?Einstein condensates</span>”, <i><span class="bib-journal">New J. Phys.</span></i>, pp. <span class="bib-pages">043030</span>, <span class="bib-year">2009</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-18"><span class="bib-index">18</span></a>] </span> <span class="bib-authors">D. Hochstuhl, S. Bauch, M. Bonitz</span>: “<span class="bib-title">Multiconfigurational time-dependent Hartree-Fock calculations for photoionization of one-dimensional Helium</span>”, <i><span class="bib-journal">Journal of Physics: Conference Series</span></i>, pp. <span class="bib-pages">012019</span>, <span class="bib-year">2010</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-19"><span class="bib-index">19</span></a>] </span> <span class="bib-authors">D. Hochstuhl, M. Bonitz</span>: “<span class="bib-title">Two-photon ionization of helium studied with the multiconfigurational time-dependent Hartree-Fock method</span>”, <i><span class="bib-journal">J. Chem. Phys.</span></i>, pp. <span class="bib-pages">084106</span>, <span class="bib-year">2011</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-20"><span class="bib-index">20</span></a>] </span> <span class="bib-authors">J. Caillat, J. Zanghellini, M. Kitzler, O. Koch, W. Kreuzer, A. Scrinzi</span>: “<span class="bib-title">Correlated multielectron systems in strong laser fields: A multiconfiguration time-dependent Hartree-Fock approach</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">012712</span>, <span class="bib-year">2005</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-21"><span class="bib-index">21</span></a>] </span> <span class="bib-authors">A. Griesmaier J. Stuhler, L. Santos</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">PRL</span></i>, pp. <span class="bib-pages">150406</span>, <span class="bib-year">2005</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-22"><span class="bib-index">22</span></a>] </span> <span class="bib-authors">J. Zanghellini, M. Kitzler, C. Fabian, T. Brabec, A. Scrinzi</span>: “<span class="bib-title">An MCTDHF Approach to Multielectron Dynamics in Laser Fields</span>”, <i><span class="bib-journal">Laser Phys.</span></i>, pp. <span class="bib-pages">1064</span>, <span class="bib-year">2003</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-23"><span class="bib-index">23</span></a>] </span> <span class="bib-authors">Jun Yan</span>: “<span class="bib-title">Harmonic Interaction Model and Its Applications in Bose-Einstein Condensation</span>”, <i><span class="bib-journal">J. Stat. Phys.</span></i>, pp. <span class="bib-pages">623</span>, <span class="bib-year">2003</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-24"><span class="bib-index">24</span></a>] </span> <span class="bib-authors">R. Kosloff</span>: “<span class="bib-title">Time-Dependent Quantum-Mechanical Mthods for Molecular-Dynamics</span>”, <i><span class="bib-journal">J. Phys. Chem.</span></i>, pp. <span class="bib-pages">2087</span>, <span class="bib-year">1988</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-25"><span class="bib-index">25</span></a>] </span> <span class="bib-authors">Leon Cohen, Chongmoon Lee</span>: “<span class="bib-title">Exact reduced density matrices for a model problem</span>”, <i><span class="bib-journal">J. Math. Phys.</span></i>, pp. <span class="bib-pages">3105</span>, <span class="bib-year">1985</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-26"><span class="bib-index">26</span></a>] </span> <span class="bib-authors">N. Q. Burdick M. Lu, B. L. Lev</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">PRL</span></i>, pp. <span class="bib-pages">190401</span>, <span class="bib-year">2011</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-27"><span class="bib-index">27</span></a>] </span> <span class="bib-authors">M. Nest, T. Klamroth, P. Saalfrank</span>: “<span class="bib-title">The multiconfiguration time-dependent Hartree-Fock method for quantum chemical calculations</span>”, <i><span class="bib-journal">J. Chem. Phys.</span></i>, pp. <span class="bib-pages">124102</span>, <span class="bib-year">2005</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-28"><span class="bib-index">28</span></a>] </span> <span class="bib-authors">H.-D. Meyer, U. Manthe, L. S. Cederbaum</span>: “<span class="bib-title">The multi-configurational time-dependent Hartree approach</span>”, <i><span class="bib-journal">Chem. Phys. Lett.</span></i>, pp. <span class="bib-pages">73—78</span>, <span class="bib-year">1990</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-29"><span class="bib-index">29</span></a>] </span> <span class="bib-authors">O. E. Alon, A. I. Streltsov, L. S. Cederbaum</span>: “<span class="bib-title">The multi-configurational time-dependent Hartree method for bosons: Many-body dynamics of bosonic systems</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">033613-1 - 4</span>, <span class="bib-year">2008</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-30"><span class="bib-index">30</span></a>] </span> <span class="bib-authors">K. Sakmann, A. I. Streltsov, O. E. Alon, L. S. Cederbaum</span>: “<span class="bib-title">Exact Quantum Dynamics of a Bosonic Josephson Junction</span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">220601-1 - 4</span>, <span class="bib-year">2009</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-31"><span class="bib-index">31</span></a>] </span> <span class="bib-authors">K. Sakmann, A. I. Streltsov, O. E. Alon, L. S. Cederbaum</span>: “<span class="bib-title">Optimal time-dependent lattice models for nonequilibrium dynamics</span>”, <i><span class="bib-journal">New J. Phys.</span></i>, pp. <span class="bib-pages">043003</span>, <span class="bib-year">2011</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-32"><span class="bib-index">32</span></a>] </span> <span class="bib-authors">Kaspar Sakmann, Alexej I. Streltsov, Ofir E. Alon, Lorenz S. Cederbaum</span>: “<span class="bib-title">Quantum dynamics of attractive versus repulsive bosonic Josephson junctions: Bose-Hubbard and full-Hamiltonian results</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">013620</span>, <span class="bib-year">2010</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-33"><span class="bib-index">33</span></a>] </span> <span class="bib-authors">F. Schreck, L. Khaykovich, K. L. Corwin, G. Ferrari, T. Bourdel, J. Cubizolles, C. Salomon</span>: “<span class="bib-title">Quasipure Bose-Einstein Condensate Immersed in a Fermi Sea</span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">080403</span>, <span class="bib-year">2001</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-34"><span class="bib-index">34</span></a>] </span> <span class="bib-authors">A. I. Streltsov, O. E. Alon, L. S. Cederbaum</span>: “<span class="bib-title">General variational many-body theory with complete self-consistency for trapped bosonic systems</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">063626</span>, <span class="bib-year">2006</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-35"><span class="bib-index">35</span></a>] </span> <span class="bib-authors">A. I. Streltsov, O. E. Alon, L. S. Cederbaum</span>: “<span class="bib-title">Role of excited states in the splitting of a trapped interacting Bose-Einstein condensate by a time-dependent barrier</span>”, <i><span class="bib-journal">Phys. Rev. Lett.</span></i>, pp. <span class="bib-pages">030402-1 - 4</span>, <span class="bib-year">2007</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-36"><span class="bib-index">36</span></a>] </span> <span class="bib-authors">A. I. Streltsov, K. Sakmann, O. E. Alon, L. S. Cederbaum</span>: “<span class="bib-title">Accurate multi-boson long-time dynamics in triple-well periodic traps</span>”, <i><span class="bib-journal">Phys. Rev. A</span></i>, pp. <span class="bib-pages">043604</span>, <span class="bib-year">2011</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-37"><span class="bib-index">37</span></a>] </span> <span class="bib-authors">T. Jun Park, J. C. Light</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">J. Chem. Phys</span></i>, pp. <span class="bib-pages">5870</span>, <span class="bib-year">1986</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-38"><span class="bib-index">38</span></a>] </span> <span class="bib-authors">Tsuyoshi Kato, Hirohiko Kono</span>: “<span class="bib-title">Time-dependent multiconfiguration theory for electronic dynamics of molecules in an intense laser field</span>”, <i><span class="bib-journal">Chem. Phys. Lett.</span></i>, pp. <span class="bib-pages">533</span>, <span class="bib-year">2004</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-39"><span class="bib-index">39</span></a>] </span> <span class="bib-authors">U. Manthe, H.-D. Meyer, L. S. Cederbaum</span>: “<span class="bib-title">Wave-packet dynamics within the multiconfiguration Hartree framework: General aspects and application to NOCl</span>”, <i><span class="bib-journal">J. Chem. Phys.</span></i>, pp. <span class="bib-pages">3199-1 - 15</span>, <span class="bib-year">1992</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-40"><span class="bib-index">40</span></a>] </span> <span class="bib-authors">O. Vendrell, M. Brill, F. Gatti, D. Lauvergnat, H.-D. Meyer</span>: “<span class="bib-title">Full dimensional (15-dimensional) quantum-dynamical simulation of the protonated water-dimer III: Mixed Jacobi-valence parametrization and benchmark results for the zero point energy, vibrationally excited states, and infrared spectrum</span>”, <i><span class="bib-journal">J. Chem. Phys.</span></i>, pp. <span class="bib-pages">234305</span>, <span class="bib-year">2009</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-41"><span class="bib-index">41</span></a>] </span> <span class="bib-authors">W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery</span>: <i><span class="bib-title">Numerical Recipes in Fortran</span></i>. <span class="bib-publisher">Cambridge University Press, Cambridge, England</span>, <span class="bib-year">1992</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-42"><span class="bib-index">42</span></a>] </span> <span class="bib-authors">W. Ketterle</span>: “<span class="bib-title">Nobel Lecture: When atoms behave as waves: Bose-Einstein condensation and the atom laser</span>”, <i><span class="bib-journal">Rev. Mod. Phys.</span></i>, pp. <span class="bib-pages">1131-1151</span>, <span class="bib-year">2002</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-43"><span class="bib-index">43</span></a>] </span> <span class="bib-authors">Inc. Wolfram Research</span>: <i><span class="bib-title">Mathematica Edition: Version 7.0</span></i>. <span class="bib-publisher">Wolfram Research, Inc., Champaign, Illinois</span>, <span class="bib-year">2008</span>.
</p>
<p class="biblio">
<span class="entry">[<a class="biblioentry" name="biblio-44"><span class="bib-index">44</span></a>] </span> <span class="bib-authors">Z. Bačić, J. C. Light</span>: “<span class="bib-title"></span>”, <i><span class="bib-journal">J. Chem. Phys</span></i>, pp. <span class="bib-pages">4594</span>, <span class="bib-year">1986</span>.
</p>

</div>
<div class="Indented">
<p><br/>
</p>

</div>
<div class="Indented">
<div class="float">
<a class="Label" name="Figure-1"> </a><div class="figure" style="max-width: 100%;">
<img class="figure" src="Fig1.png" alt="figure Fig1.png" style="max-width: 1650px; max-height: 1275px;"/>
 <div class="caption">
Figure 1  (color online). Numerical convergence of the self-consistent MCTDHB and fixed-orbital full configuration interaction (FCI) methods for the ground state energy of the harmonic interaction model (HIM). Systems with N=2,10,50,100 and 1000 bosons are considered, the strengths of the interparticle interactions <span class="formula"><i>K</i><sub>0</sub></span> have been chosen to keep <span class="formula">Λ = <i>K</i><sub>0</sub>(<i>N</i> − 1) = 0.5</span> constant. We plot the relative differences between the total energies computed using the MCTDHB (filled symbols) and FCI (open symbols) many-body methods and respective exact energies in percents, <span class="formula">100⋅(<i>E</i><sub><span class="mathrm">MB</span></sub> − <i>E</i><sub><span class="mathrm">exact</span></sub>) ⁄ <i>E</i><sub><span class="mathrm">exact</span></sub></span>, for different orbital number <span class="formula"><i>M</i></span>. For a given <span class="formula"><i>M</i></span> both many-body methods span the same Fock space, i.e., the respective secular matrices to be diagonalized are of the same size. The advantage of the appropriate, i.e., self-consistent, choice of the one-particle basis functions is evident — the self-consistent MCTDHB method converges much faster than the fixed-orbital FCI one. Note the logarithmic scale and number of decades spaned. All quantaties shown are dimensionless.
</div>
<div class="PlainVisible">
<a class="Label" name="HIM_MCTDHB"> </a> 
</div>

</div>

</div>

</div>
<div class="Indented">
<div class="float">
<a class="Label" name="Figure-2"> </a><div class="figure">
<img class="embedded" src="Fig2a.png" alt="figure Fig2a.png" style="max-width: 1650px; max-height: 1275px;"/>
 <img class="embedded" src="Fig2b.png" alt="figure Fig2b.png" style="max-width: 1650px; max-height: 1275px;"/>
 <div class="caption">
Figure 2  (color online). A sudden change (quench) of the interparticle interaction leads to &ldquo;breathing&rdquo; dynamics of the system. We study the HIM system with <span class="formula"><i>N</i> = 2</span> bosons where the interparticle interaction strength is quenched from zero to <span class="formula"><i>K</i><sub>0</sub> = 0.5</span>. The evolution of the one-particle density at the origin <span class="formula"><i>ρ</i>(<i>x</i> = <i>x</i>’ ≡ 0)</span> is plotted as a function of time. The exact dynamics reveals oscillations with the main breathing frequency <span class="formula"><i>ω</i><sub><i>breath</i></sub> = 2<span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><sup>2</sup> + 4<i>K</i><sub>0</sub></span><span class="ignored">)</span></span></span> augmented by overtones <span class="formula">2<i>ω</i><sub><i>breath</i></sub></span>, <span class="formula">3<i>ω</i><sub><i>breath</i></sub>, …</span> (see text for more details). The solid (black) line depicts the guiding <span class="formula"> ~ cos(<i>w</i><sub><i>breath</i></sub><i>t</i>)</span> function. The numerical MCTDHB and FCI results are contrasted and compared with the exact ones, plotted by a bold (red) line. The left panel depicts the density oscillations at short times. At the FCI level accurate description of the dynamics is achieved by using at least eight fixed-shape orbitals. To gain similar accuracy within the MCTDHB method one needs only three time-adaptive orbitals. Four time-adaptive orbitals [MCTDHB(<span class="formula"><i>M</i> = 4</span>)] provide a numerically exact description. The right panel shows the density oscillations at longer times. To describe the dynamics in this case a larger Fock space (more orbitals) is required. The numerically exact description is obtained by using six time-adaptive [MCTDHB(6)] or twelve fixed-shape orbitals [FCI(12)]. See text for further discussion. All quantities shown are dimensionless.
</div>
<div class="PlainVisible">
<a class="Label" name="HIM_QUENCH_N2"> </a> 
</div>

</div>

</div>

</div>
<div class="Indented">
<div class="float">
<a class="Label" name="Figure-3"> </a><div class="figure" style="max-width: 100%;">
<img class="figure" src="Fig3.png" alt="figure Fig3.png" style="max-width: 1650px; max-height: 1275px;"/>
 <div class="caption">
Figure 3 (color online). Breathing dynamics of the HIM system with <span class="formula"><i>N</i> = 10</span> for the same interaction quench scenario as in Fig. <a class="Reference" href="#HIM_QUENCH_N2">2↑</a>. The evolution of the one-particle density at the origin <span class="formula"><i>ρ</i>(<i>x</i> = <i>x</i>’ ≡ 0)</span> is plotted, notice different scales for the short and long times. The density oscillation is formed by the main breathing frequency <span class="formula"><i>ω</i><sub><i>breath</i></sub> = 2<span class="sqrt"><span class="radical">√</span><span class="ignored">(</span><span class="root"><i>ω</i><sup>2</sup> + 2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span></span> with strong contributions of the overtones <span class="formula">2<i>ω</i><sub><i>breath</i></sub>, 3<i>ω</i><sub><i>breath</i></sub>, ...</span> (see text for more details). The solid (black) line depicts the guiding <span class="formula"> ~ cos(<i>w</i><sub><i>breath</i></sub><i>t</i>)</span> function. The MCTDHB(<span class="formula"><i>M</i> = 8</span>) method with eight time-adaptive orbitals provides very accurate description of the breathing dynamics for the short and long times. The MCTDHB results for <span class="formula"><i>M</i> = 9</span> and <span class="formula"><i>M</i> = 10</span> are identical, indicating that the exact description has been numerically reached. The FCI(<span class="formula"><i>M</i> = 16</span>) results plotted by triangles start to deviate from the exact solution already for short times. The FCI method with sixteen fixed-shape orbitals provides a reasonable description of the dynamics for a very short time only, i.e., it is incapable to describe more than a half of the first breathing cycle. The exact results could not be obtained in this model from the analytical solution Eq. (<a class="Reference" href="#TDEXPAND">7↑</a>) because it is more difficult to perform the needed 10-dimensional integrations than to solve the problem numerically exactly by MCTDHB. All quantities shown are dimensionless.
</div>
<div class="PlainVisible">
<a class="Label" name="HIM_QUENCH_N10"> </a> 
</div>

</div>

</div>

</div>
<div class="Indented">
<div class="float">
<a class="Label" name="Figure-4"> </a><div class="figure">
<img class="embedded" src="Fig4_f1.png" alt="figure Fig4_f1.png" style="max-width: 1650px; max-height: 1275px;"/>
 <div class="caption">
Figure 4  (color online). The HIM model with time-dependent trap <span class="formula"><i>V</i>(<i>x</i>) = <i>ω</i><sub>0</sub><span class="symbol">[</span>1 + <i>f</i>(<i>t</i>)<span class="symbol">]</span><i>x</i><sup>2</sup></span> and time-dependent interparticle interaction <span class="formula"><i>W</i>(<i>x</i><sub><i>i</i></sub> − <i>x</i><sub><i>j</i></sub>) = <i>K</i><sub>0</sub><span class="array"><span class="arrayrow"><span class="bracket align-left">⎡</span></span><span class="arrayrow"><span class="bracket align-left">⎣</span></span></span>1 − <span class="fraction"><span class="ignored">(</span><span class="numerator"><i>ω</i><span class="scripts"><sup class="script">2</sup><sub class="script">0</sub></span></span><span class="ignored">)/(</span><span class="denominator">2<i>NK</i><sub>0</sub></span><span class="ignored">)</span></span><i>f</i>(<i>t</i>)<span class="array"><span class="arrayrow"><span class="bracket align-right">⎤</span></span><span class="arrayrow"><span class="bracket align-right">⎦</span></span></span>(<i>x</i><sub><i>i</i></sub> − <i>x</i><sub><i>j</i></sub>)<sup>2</sup></span> permits exact solution. The exact expectation value of the total Hamiltonian of the system Eq. (<a class="Reference" href="#TDEXPVALE">13↑</a>) reads <span class="formula">⟨Ψ(<i>t</i>)∣<i>Ĥ</i>(<i>t</i>)∣Ψ(<i>t</i>)⟩ = <i>ϵ</i>(<i>t</i>) + <i>const</i>.</span>, with the time-independent constant equals to <span class="formula"><span class="fraction"><span class="ignored">(</span><span class="numerator"><i>D</i></span><span class="ignored">)/(</span><span class="denominator">2</span><span class="ignored">)</span></span>(<i>N</i> − 1)<i>δ</i><sub><i>N</i></sub></span> and <span class="formula"><i>D</i> = 1</span>. The driven function <span class="formula"><i>f</i>(<i>t</i>) = <i>f</i><sub>1</sub>(<i>t</i>)</span> and the time-dependent part of the energy <span class="formula"><i>ϵ</i>(<i>t</i>) = <i>ϵ</i><sub>1</sub>(<i>t</i>)</span> for N=10 bosons with <span class="formula"><i>K</i><sub>0</sub> = 0.5</span> are plotted. The convergence of <span class="formula"><i>ϵ</i><sub>1</sub>(<i>t</i>)</span> when increasing the number of the time-adaptive orbitals <span class="formula"><i>M</i></span> is depicted. The Gross-Pitaevskii results [GP <span class="formula"> ≡ </span>MCTDHB(1)], plotted by a bold solid line, are inaccurate even for very short time. The MCTDHB(3) provides excellent description up to <span class="formula"><i>t</i> ≈ 5</span>, the MCTDHB(5) works well till <span class="formula"><i>t</i> ≈ 15</span>, the MCTDHB(6) till <span class="formula"><i>t</i> ≈ 30</span>, and the MCTDHB(7) results coincide with the exact solution at all the times depicted. See text for discussion. All quantities shown are dimensionless.
</div>
<div class="PlainVisible">
<a class="Label" name="HIMTD_ft1"> </a> 
</div>

</div>

</div>

</div>
<div class="Indented">
<div class="float">
<a class="Label" name="Figure-5"> </a><div class="figure">
<img class="embedded" src="Fig4_f2.png" alt="figure Fig4_f2.png" style="max-width: 1650px; max-height: 1275px;"/>
<div class="caption">
Figure 5 (color online). The modified HIM model with time-dependent trap and interparticle interaction driven by a complicated function <span class="formula"><i>f</i><sub>2</sub>(<i>t</i>)</span> [Eq. (<a class="Reference" href="#f1_f2">14↑</a>)]. The function is depicted in the upper panel. The time-dependent contribution <span class="formula"><i>ϵ</i><sub>2</sub>(<i>t</i>)</span> to the total energies is computed at several different levels of the MCTDHB(<span class="formula"><i>M</i></span>) theory for <span class="formula"><i>N</i> = 10</span> [<span class="formula"><i>M</i> = 6, 7</span>] and <span class="formula"><i>N</i> = 50</span> [<span class="formula"><i>M</i> = 5, 6</span>] bosons. The strength of the interparticle interaction is <span class="formula"><i>K</i><sub>0</sub> = 0.5</span>. The considered time-dependency of the one- and two-body interaction potentials guarantees that the exact <span class="formula"><i>ϵ</i><sub>2</sub>(<i>t</i>)</span>, plotted by open red circles, is the same for both systems. The MCTDHB(5) for <span class="formula"><i>N</i> = 50</span> and MCTDHB(6) for <span class="formula"><i>N</i> = 10</span> provide converged description of the dynamics till <span class="formula"><i>t</i> ≈ 25</span>, for longer times more orbitals are needed for &ldquo;absolute&rdquo; convergence. The corresponding Gross-Pitaevskii results, marked by arrows, are semi-qualitative for very short initial times only, till <span class="formula"><i>t</i> ≈ 1</span>. See text for discussion. All quantities shown are dimensionless.
</div>
<a class="Label" name="HIMTD_ft2"> </a> 
</div>

</div>

</div>

<hr class="footer"/>
<div class="footer" id="generated-by">
Document generated by <a href="http://elyxer.nongnu.org/">eLyXer 1.2.2 (2011-06-12)</a> on <span class="create-date">2015-04-01T13:13:28.999192</span>
</div>
</div>
</body>
</html>
