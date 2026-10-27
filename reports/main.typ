/////////////
// Imports //
/////////////
#import "@preview/adaptable-pset:0.2.0": *
#import "@preview/physica:0.9.8": *
#import "@preview/unify:0.8.1": *
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#show: codly-init.with()
#codly(languages: codly-languages)

/////////////////
// Maths Setup //
/////////////////

// upright vectors
#let vectorboldupright(a) = vb($upright(#a)$)
#let vbu = vectorboldupright
#let vectorunitupright(a) = vu($upright(#a)$)
#let vuu = vectorunitupright
#let vectorarrowupright(a) = va($upright(#a)$)
#let vau = vectorarrowupright

// automatically use square brackets for vectors and matrices
#set math.vec(delim: "[")
#set math.mat(delim: "[")
#let vecrowOld = vecrow
#let vecrow = vecrowOld.with(delim: "[")

// custom units
#add-unit("decibelIsotropic", "dBi", "upright(\"dBi\")", space: true)
#add-unit("decibelWatt", "dBW", "upright(\"dBW\")", space: true)

////////////////////
// Document Setup //
////////////////////

// assignment info
#show: homework.with(
    title: "HW03",
    author: "Vai Srivastava",
    collaborators: [],
    course-id: "ENAE 641: Linear System Dynamics",
    instructor: "Dr. Robert Sanner",
    semester: "Fall 2026",
    due-time: datetime(
        // due date
        year: 2026,
        month: 10,
        day: 27,

        // due time
        hour: 23,
        minute: 59,
        second: 00,
    ).display(
        "[month repr:long] [day padding:zero], [year] at [hour repr:24]:[minute padding:zero]:[second padding:zero]",
    ),

    // (defaults to A4)
    paper-size: "us-letter",
)

// document settings
#set text(font: "New Computer Modern", size: 10pt)
#set enum(numbering: "a)")

// problem headings
#let probOld = prob
#let prob = prob.with(color: black)

////////////////////////////
// The Assignment Itself: //
// Problems and Solutions //
////////////////////////////

#prob[
    The file `F16.mat` on ELMS contains a standard state-space perturbation model for the longitudinal dynamics of an F16 aircraft in level flight at #qty(1e5, "ft") moving with speed #qty(500, "ft/s"). The four states of this model are: speed, angle of attack, pitch angle, and pitch rate, in that order.

    + Determine the real block modal form for the $vb(A)$ matrix, and show that the system is stable.
    <hwk:p01a>

    + Specify a physical (i.e. real-valued) initial condition which will uniquely excite the phugoid (slowly oscillating) mode. Repeat for the short period (fast oscillating) mode. Interpret each in terms of the physical variables.
    <hwk:p01b>

    + The two inputs to this model are throttle and elevator deviations from trim, respectively. Which mode is most strongly excited by the throttle? Which mode is most strongly excited by the elevator?

        #smallcaps[Hint:] Look at the corresponding modal $vb(B)$ matrix for the transform in #link(<hwk:p01a>)[Part A]
    <hwk:p01c>
] <hwk:p01>

+ Answer
<hwk:s01a>

+ Answer
<hwk:s01b>

+ Answer
<hwk:s01c>

#pagebreak(weak: true)

#prob[
    Consider a time-periodic system from #link("https://github.com/vaisriv/enae641-hw02")[HW02]:
    $
        vb(A)(t) = mat(
            -1 + 1.5 cos^2(t), 1 - 1.5 sin(t) cos(t);
            -1 - 1.5 sin(t) cos(t), 1 + 1.5 sin^2(t)
        )
    $
    This system is time-periodic, hence it admits a Floquet decomposition. Determine the components $(vb(P)(t), vb(Q), vb(R))$ of this decomposition.
] <hwk:p02>

Answer
<hwk:s02>

#pagebreak(weak: true)

#prob[
    For the linear time-varying system
    $
        dot(vbu(x))(t) = vb(A)(t) vbu(x)(t) quad vbu(x)(t_0) = vbu(x)_0
    $
    if
    $
        vb(A)(t) = mat(
            0, t sin(t) cos(sqrt(t));
            -t sin(t) cos(sqrt(t)), 0;
        )
    $
    and $vbu(x)_0^TT = vecrow(-1, 2)$, find $||vbu(x)(t)||_2$ for all $t >= t_0$. What does this tell you about the stability properties of this system (unstable, asymptotically stable, or stable)?
] <hwk:p03>

Answer
<hwk:s03>

#pagebreak(weak: true)

#prob[
    + For the F16 model in #link(<hwk:p01>)[Problem 1], find a constant, positive-definite matrix $vb(P)$ so that
        $
            dv((vbu(x)^TT(t) P vbu(x)(t)), t, style: "large") = - ||vbu(x)(t)||^2
        $
        (Use `MATLAB`'s ```matlab lyap()``` function, but read its documentation carefully!) Verify this $vb(P)$ is positive definite.
    <hwk:p04a>

    + Use this $(vb(P), vb(Q))$ pair to bound the exponential convergence of $||vbu(x)(t)||$ from an arbitrary $vbu(x)_0$.
    <hwk:p04b>
] <hwk:p04>

+ Answer
<hwk:s04a>

+ Answer
<hwk:s04b>

#pagebreak(weak: true)

#prob[
    + For the F16 model in #link(<hwk:p01>)[Problem 1], what are the eigenvalues of $vb(A) + vb(A)^TT$? Do these eigenvalues allow you to conclude the kind of exponential stability result that you obitained in #link(<hwk:p01>)[Problem 1] and #link(<hwk:p04>)[Problem 4]?
    <hwk:p05a>

    + Is there a contradiction between the results of #link(<hwk:p05a>)[Part A] and that of #link(<hwk:p01>)[Problem 1] and #link(<hwk:p04>)[Problem 4]? Explain carefully in light of Lyapunov theory.

        #smallcaps[Hint:] What $vb(P)$ are you implicitly using in the test of #link(<hwk:p05a>)[Part A]? Should you expect an arbitrary choice such as this to satisfy Lyapunov's requirement?
    <hwk:p05b>
] <hwk:p05>

+ Answer
<hwk:s05a>

+ Answer
<hwk:s05b>

#pagebreak(weak: true)

#prob[
    For the homogenous dynamics of a linear system, suppose that the symmetric matrices $vb(P)(t), vb(Q)(t)$ satisfy the Lyapunov equation, with $vb(P)(t)$ bounded and uniformly positive definite, and $vb(Q)(t)$ bounded. Let $nu(t)$ be the minimal eigenvalue of $vb(Q)(t)$, i.e.
    $
        nu(t) = min_k lambda_k (vb(Q)(t))
    $
    Determine the "strongest" property (in order: inconclusive, stable, asymptotically stable, or exponentially stable) you can conclude for the following cases below. Show your reasoning carefully, and note that some or all of these cases might not fit neatly into the categories specifically considered in class.

    + $nu(t) = 1/t$
    <hwk:p06a>

    + $nu(t) = -c e^(-t)$ with $c > 0$
    <hwk:p06b>

    + $nu(t) = 1/t^2$
    <hwk:p06c>

    + $nu(t) = 2 - 3 sin(t)$
    <hwk:p06d>
] <hwk:p06>

+ Answer
<hwk:s06a>

+ Answer
<hwk:s06b>

+ Answer
<hwk:s06c>

+ Answer
<hwk:s06d>

#pagebreak(weak: true)

#prob[
    Generalize your observations from #link(<hwk:p06>)[Problem 6]. What general conditions on $nu(t)$ would allow you to conclude isL stability? Asymptotic stability? Exponential stability?

    #smallcaps[Note:] This question is asking you to go beyond the comparatively straightforward results discussed in lecture; for example to demonstrate formally that stability can sometimes be concluded if $nu(t) < 0$ for some or even all $t >= t_0$, etc.
] <hwk:p07>

Answer
<hwk:s07>

#pagebreak(weak: true)

== Code

#codly(header: [./src/index.m])
#raw(read("../src/index.m"), block: true, lang: "matlab") <code:index.m>
