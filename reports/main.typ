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
    The file `F16.mat` on ELMS contains a standard state-space perturbation model for the longitudinal dynamics of an F16 aircraft in level flight at #qty(1e4, "ft") moving with speed #qty(500, "ft/s"). The four states of this model are: speed, angle of attack, pitch angle, and pitch rate, in that order.

    + Determine the real block modal form for the $vb(A)$ matrix, and show that the system is stable.
    <hwk:p01a>

    + Specify a physical (i.e. real-valued) initial condition which will uniquely excite the phugoid (slowly oscillating) mode. Repeat for the short period (faster oscillating) mode. Interpret each in terms of the physical variables.
    <hwk:p01b>

    + The two inputs to this model are throttle and elevator deviations from trim, respectively. Which mode is most strongly excited by the throttle? Which mode is most strongly excited by the elevator?

        #smallcaps[Hint:] Look at the corresponding modal $vb(B)$ matrix for the transform in #link(<hwk:p01a>)[Part A]
    <hwk:p01c>
] <hwk:p01>

+
    The eigenvalues of $vb(A)$ form two complex-conjugate pairs:
    $
        lambda_(1,2) = sigma_1 plus.minus j omega_1 = -0.01425 plus.minus 0.04434 j, quad
        lambda_(3,4) = sigma_3 plus.minus j omega_3 = -1.695 plus.minus 2.058 j
    $
    The slow pair is the phugoid ($omega_n = #qty(0.04658, "rad/s")$, $zeta = 0.306$, period $2 pi\/omega_1 approx #qty(141.7, "s")$). The fast pair is the short period ($omega_n = #qty(2.666, "rad/s")$, $zeta = 0.636$, period $2 pi\/omega_3 approx #qty(3.05, "s")$).

    Let $vbu(v)_1$ and $vbu(v)_3$ be the eigenvectors for $lambda_1 = sigma_1 + j omega_1$ and $lambda_3 = sigma_3 + j omega_3$. Taking real and imaginary parts of $vb(A) vbu(v)_1 = lambda_1 vbu(v)_1$ gives
    $
        vb(A) mat(Re{vbu(v)_1}, Im{vbu(v)_1}) = mat(Re{vbu(v)_1}, Im{vbu(v)_1}) mat(sigma_1, omega_1; -omega_1, sigma_1)
    $
    and likewise for $vbu(v)_3$. So the real transformation $vbu(x) = vb(P) vbu(z)$ with
    $
        vb(P) = mat(Re{vbu(v)_1}, Im{vbu(v)_1}, Re{vbu(v)_3}, Im{vbu(v)_3})
    $
    puts $vb(A)_p = vb(P)^(-1) vb(A) vb(P)$ in block modal form, with the phugoid in $(z_1, z_2)$ and the short period in $(z_3, z_4)$. Each eigenvector is scaled so that its largest entry is $1$: $vbu(v)_1$ by its speed entry, and $vbu(v)_3$ by its pitch-rate entry. `MATLAB` then gives
    $
        vb(P) = mat(
            1, 0, 0.01390, -0.01653;
            -0.008132, -0.0001973, -0.06345, -0.4306;
            -0.2367, -0.6968, -0.1014, -0.1231;
            0.08060, -0.001325, 1, 0;
        )
    $
    $
        vb(A)_p = vb(P)^(-1) vb(A) vb(P) = mat(
            -0.01425, 0.04434, 0, 0;
            -0.04434, -0.01425, 0, 0;
            0, 0, -1.695, 2.058;
            0, 0, -2.058, -1.695;
        ) quad qed
    $

    #emph[Stability:] $vb(A)_p$ is block diagonal, so with the rotation matrix $vb(Theta)(theta) = mat(cos(theta), sin(theta); -sin(theta), cos(theta))$,
    $
        e^(vb(A)_p t) = mat(e^(sigma_1 t) vb(Theta)(omega_1 t), vb(0); vb(0), e^(sigma_3 t) vb(Theta)(omega_3 t)), quad
        vbu(x)(t) = vb(P) e^(vb(A)_p t) vb(P)^(-1) vbu(x)_0
    $
    Each block is a rotation (norm $1$) scaled by $e^(sigma_k t)$, and $sigma_3 < sigma_1 < 0$, so $norm(e^(vb(A)_p t))_(i,2) = e^(sigma_1 t)$ and
    $
        norm(vbu(x)(t))_2 <= norm(vb(P))_(i,2) norm(vb(P)^(-1))_(i,2) e^(sigma_1 t) norm(vbu(x)_0)_2 -> 0
    $
    for every $vbu(x)_0$. Both modes decay exponentially, i.e. all poles lie in the open left half-plane:
    $
        max_k Re(lambda_k) = sigma_1 = -0.01425 < 0 quad => quad "asymptotically (exponentially) stable" quad qed
    $
    The `MATLAB` calculation:
    #raw(read("../outputs/text/s01a.txt"), block: true)
<hwk:s01a>

+
    An initial state in a mode's plane excites only that mode. For $vbu(x)_0 = a Re{vbu(v)_1} + b Im{vbu(v)_1}$ with real $a$ and $b$, $vbu(z)_0 = vb(P)^(-1) vbu(x)_0 = vecrow(a, b, 0, 0)^TT$. Since $e^(vb(A)_p t)$ is block diagonal, $z_3(t) = z_4(t) = 0$ for all $t$, and
    $
        vbu(x)(t) = e^(sigma_1 t) [a cos(omega_1 t) + b sin(omega_1 t)] Re{vbu(v)_1} + e^(sigma_1 t) [b cos(omega_1 t) - a sin(omega_1 t)] Im{vbu(v)_1}
    $
    contains only the phugoid and stays in the plane $op("span"){Re{vbu(v)_1}, Im{vbu(v)_1}}$. Any real $a$ and $b$ that are not both zero work. Since the model is linear, the initial state can be scaled down to stay a small perturbation from trim. Taking $a = 1$ and $b = 0$:
    $
        vbu(x)_0^"ph" = Re{vbu(v)_1} = vec(1, -0.008132, -0.2367, 0.08060) quad qed
    $
    Likewise, $vbu(x)_0 = Re{vbu(v)_3}$ gives $vbu(z)_0 = vecrow(0, 0, 1, 0)^TT$ and excites only the short period:
    $
        vbu(x)_0^"sp" = Re{vbu(v)_3} = vec(0.01390, -0.06345, -0.1014, 1) quad qed
    $

    #emph[Interpretation:] Write each eigenvector in polar form, $v_i = r_i e^(j phi_i)$. While a mode is excited, state $x_i$ oscillates with relative amplitude $r_i$ and relative phase $phi_i$:
    #align(center, table(
        columns: 5,
        align: (left, right, right, right, right),
        stroke: none,
        table.hline(),
        table.header(
            [state],
            [phugoid $r_i$],
            [phugoid $phi_i$],
            [short period $r_i$],
            [short period $phi_i$],
        ),
        table.hline(stroke: 0.5pt),
        [speed],
        [$1$],
        [#qty(0, "degree")],
        [$0.0216$],
        [#qty(-49.9, "degree")],
        [angle of attack],
        [$0.00813$],
        [#qty(-178.6, "degree")],
        [$0.435$],
        [#qty(-98.4, "degree")],
        [pitch angle],
        [$0.736$],
        [#qty(-108.8, "degree")],
        [$0.160$],
        [#qty(-129.5, "degree")],
        [pitch rate],
        [$0.0806$],
        [#qty(-0.9, "degree")],
        [$1$],
        [#qty(0, "degree")],
        table.hline(),
    ))
    The output matrix $vb(C) = mat(500, 0, 0, 0; 0, -500, 500, 0)$ suggests that the speed state is normalized by $V_0 = #qty(500, "ft/s")$. Then $vbu(y)$ holds the speed change and the climb rate $V_0 (theta - alpha)$, both in #unit("ft/s"), and $r_i$ compares a fractional speed change with $alpha$ and $theta$ in radians. The pitch-rate state has its own scaling: $dot(theta) = 0.4252 thin q$ from the third row of $vb(A)$. So its $r_i$ is in those units, and measured as $dot(theta)$ the short period's pitch-rate amplitude is about $0.43$, comparable to that of $alpha$.
    - #emph[Phugoid:] speed and pitch angle carry the motion. Pitch angle lags speed by about #qty(109, "degree"), while the angle of attack barely moves ($r = 0.008$). With $alpha$ nearly constant, the flight path angle $theta - alpha approx theta$, so the aircraft climbs as it slows down and dives as it speeds up. It trades kinetic for potential energy over a slow (about #qty(142, "s")) oscillation at nearly constant angle of attack.
    - #emph[Short period:] pitch rate and angle of attack carry the motion. Angle of attack lags pitch rate by about #qty(98, "degree"), and speed barely changes ($r = 0.02$). This is a fast (about #qty(3, "s")), well-damped ($zeta = 0.64$) pitching of the airframe about its center of mass at nearly constant speed.

    The two initial states read the same way. $vbu(x)_0^"ph"$ starts the aircraft above trim speed with the nose down ($theta = -0.237$ per unit of fractional speed), almost no change in $alpha$, and a small pitch rate. $vbu(x)_0^"sp"$ is essentially a pitch-rate kick, with small negative offsets in $alpha$ and $theta$ and almost no change in speed.

    Simulating both initial states confirms that the other mode's modal coordinates stay at zero (to about $10^(-16)$):
    #raw(read("../outputs/text/s01b.txt"), block: true)
    #figure(
        image("../outputs/figures/s01b.png", width: 100%),
        caption: [Homogeneous response (top) and modal coordinates $vbu(z) = vb(P)^(-1) vbu(x)$ (bottom) for the phugoid initial state $Re{vbu(v)_1}$ (left) and the short-period initial state $Re{vbu(v)_3}$ (right).],
    )
<hwk:s01b>

+
    In modal coordinates $dot(vbu(z)) = vb(A)_p vbu(z) + vb(B)_p vbu(u)$ with $vb(B)_p = vb(P)^(-1) vb(B)$. So entry $(k, j)$ of $vb(B)_p$ is how strongly input $u_j$ drives $z_k$. Rows 1–2 drive the phugoid and rows 3–4 the short period; column 1 is throttle and column 2 is elevator:
    $
        vb(B)_p = vb(P)^(-1) vb(B) = mat(
            0.05535, 0.006327;
            -0.01890, 0.03786;
            -0.004487, -0.3446;
            0.004247, 0.05742;
        ) quad qed
    $
    Comparing the size of each mode's pair of rows, for each input:
    $
        "throttle:" quad & norm(vecrow(0.05535, -0.01890)) = 0.0585 quad "vs" quad norm(vecrow(-0.004487, 0.004247)) = 0.00618 quad (approx 9.5 times) \
        "elevator:" quad & norm(vecrow(0.006327, 0.03786)) = 0.0384 quad "vs" quad norm(vecrow(-0.3446, 0.05742)) = 0.349 quad (approx 9.1 times)
    $
    with the phugoid block listed first in each line. The size of $vb(B)_p$ depends on how the eigenvectors are scaled: rescaling $vbu(v)_k$ rescales mode $k$'s rows inversely. A comparison that does not depend on that scaling is the part of each input direction $vbu(b)_j$ lying in mode $k$'s plane, $vb(P)_k vb(B)_(p,k) vbu(e)_j$, where $vb(P)_k$ holds mode $k$'s two columns of $vb(P)$ and $vb(B)_(p,k)$ its two rows of $vb(B)_p$. Its norm is $0.0555$ vs $0.00475$ for the throttle ($approx 11.7$ times) and $0.0286$ vs $0.346$ for the elevator ($approx 12.1$ times). Both measures agree:
    $
        "throttle" => "phugoid", quad "elevator" => "short period" quad qed
    $
    This matches the physics, and the columns of $vb(B)$ show it directly. The throttle changes thrust, which acts mainly on speed ($vb(B)_(:,1) approx 0.055 vbu(e)_1$), and speed is the dominant component of the phugoid. The elevator produces a pitching moment, which acts mainly on pitch rate ($vb(B)_(:,2) approx -0.344 vbu(e)_4$), and pitch rate and angle of attack are what move in the short period.
    #raw(read("../outputs/text/s01c.txt"), block: true)
<hwk:s01c>

#pagebreak(weak: true)

#prob[
    Consider a time-periodic system from #link("https://github.com/vaisriv/enae641-hw02")[HW02]:
    $
        vb(A)(t) = mat(
            -1 + 1.5 cos^2(t), 1 - 1.5 sin(t) cos(t);
            -1 - 1.5 sin(t) cos(t), -1 + 1.5 sin^2(t)
        )
    $
    This system is time-periodic, hence it admits a Floquet decomposition. Determine the components $(vb(P)(t), vb(Q), vb(R))$ of this decomposition.
] <hwk:p02>

Every entry of $vb(A)(t)$ is built from $cos^2(t)$, $sin(t) cos(t)$, and $sin^2(t)$, all of which have period $pi$. So $vb(A)(t + pi) = vb(A)(t)$, and the system is LTP with $T = pi$. Take $t_0 = 0$.

#emph[State transition matrix:] Define the rotation matrix $vb(Theta)(t)$ and the constant matrices $vb(Lambda)$ and $vb(S)$:
$
    vb(Theta)(t) = mat(cos(t), sin(t); -sin(t), cos(t)), quad
    vb(Lambda) = mat(1/2, 0; 0, -1), quad
    vb(S) = mat(0, 1; -1, 0)
$
These satisfy $dot(vb(Theta)) = vb(S) vb(Theta) = vb(Theta) vb(S)$ and $vb(Theta)(t + pi) = -vb(Theta)(t)$. Multiplying out (with $s = sin(t)$, $c = cos(t)$, and $c^2 + s^2 = 1$) shows that $vb(A)(t)$ is a rotated copy of $vb(Lambda)$, plus $vb(S)$:
$
    vb(Theta) vb(Lambda) vb(Theta)^TT + vb(S)
    = mat(1/2 c^2 - s^2, -3/2 s c; -3/2 s c, 1/2 s^2 - c^2) + mat(0, 1; -1, 0)
    = mat(-1 + 1.5 c^2, 1 - 1.5 s c; -1 - 1.5 s c, -1 + 1.5 s^2)
    = vb(A)(t)
$
Substituting $vbu(x) = vb(Theta)(t) vbu(w)$ into $dot(vbu(x)) = vb(A)(t) vbu(x)$ then gives a time-invariant system:
$
    dot(vbu(w)) = vb(Theta)^TT [vb(A) vb(Theta) - dot(vb(Theta))] vbu(w)
    = vb(Theta)^TT [vb(Theta) vb(Lambda) + vb(S) vb(Theta) - vb(S) vb(Theta)] vbu(w)
    = vb(Lambda) vbu(w)
$
so the STM is
$
    vb(Phi)(t, tau) = vb(Theta)(t) e^(vb(Lambda) (t - tau)) vb(Theta)^TT (tau)
$
which agrees with the STM found in #link("https://github.com/vaisriv/enae641-hw02")[HW02].

#emph[Monodromy matrix:] Since $vb(Theta)(pi) = -vb(I)$ and $vb(Theta)(0) = vb(I)$,
$
    vb(Q) = vb(Phi)(t_0 + T, t_0) = vb(Phi)(pi, 0) = -e^(vb(Lambda) pi) = mat(-e^(pi\/2), 0; 0, -e^(-pi)) approx mat(-4.810, 0; 0, -0.0432)
$

#emph[Constant matrix $vb(R)$:] $vb(R)$ is defined by $e^(vb(R) T) = vb(Q)$, so $vb(R) = 1/T log vb(Q)$. No real $vb(R)$ works here. The eigenvalues of $vb(Q)$ are $e^(pi r_k)$, where $r_k$ are the eigenvalues of $vb(R)$. Making $e^(pi r_k)$ negative requires $Im(r_k)$ to be an odd integer, so $r_k$ is not real. A real $vb(R)$ would then have the eigenvalue pair $r_k$, $macron(r)_k$, whose exponentials are equal. That would make the eigenvalues of $vb(Q)$ equal, but $-e^(pi\/2) != -e^(-pi)$. The principal matrix logarithm (`MATLAB`'s ```matlab logm()```), with $log(-a) = ln(a) + j pi$, gives
$
    vb(R) = 1/pi log vb(Q) = 1/pi mat(pi\/2 + j pi, 0; 0, -pi + j pi) = mat(1/2 + j, 0; 0, -1 + j) = vb(Lambda) + j vb(I)
$

#emph[Periodic part $vb(P)(t)$:] Since $vb(Lambda)$ and $j vb(I)$ commute, $e^(-vb(R) t) = e^(-vb(Lambda) t) e^(-j t)$, and
$
    vb(P)(t) = vb(Phi)(t, t_0) e^(-vb(R) t) = vb(Theta)(t) e^(vb(Lambda) t) e^(-vb(Lambda) t) e^(-j t) = e^(-j t) vb(Theta)(t), quad
    vb(P)^(-1)(t) = e^(j t) vb(Theta)^TT (t)
$
This $vb(P)(t)$ has the properties the decomposition requires:
- It is periodic: $vb(P)(t + pi) = e^(-j t) e^(-j pi) vb(Theta)(t + pi) = e^(-j t) (-1)(-vb(Theta)(t)) = vb(P)(t)$.
- It is a Lyapunov transformation: $abs(det vb(P)(t)) = abs(e^(-2 j t)) det vb(Theta)(t) = 1$ and $abs(det vb(P)^(-1)(t)) = 1$ for all $t$, so both are bounded away from zero.
- It reproduces the STM: $vb(P)(t) e^(vb(R)(t - tau)) vb(P)^(-1)(tau) = e^(-j t) vb(Theta)(t) e^(vb(Lambda)(t - tau)) e^(j(t - tau)) e^(j tau) vb(Theta)^TT (tau) = vb(Phi)(t, tau)$.

So the Floquet decomposition $vb(Phi)(t, tau) = vb(P)(t) e^(vb(R)(t - tau)) vb(P)^(-1)(tau)$, with period $T = pi$, has
$
    vb(P)(t) = e^(-j t) mat(cos(t), sin(t); -sin(t), cos(t)), quad
    vb(Q) = mat(-e^(pi\/2), 0; 0, -e^(-pi)), quad
    vb(R) = mat(1/2 + j, 0; 0, -1 + j) quad qed
$

#emph[A real-valued alternative:] $vb(A)(t)$ is also $2 pi$-periodic. Over $T = 2 pi$ the monodromy matrix has positive eigenvalues, so it has a real logarithm and every component is real:
$
    vb(Q) = vb(Phi)(2 pi, 0) = e^(2 pi vb(Lambda)) = mat(e^pi, 0; 0, e^(-2 pi)), quad
    vb(R) = 1/(2 pi) log vb(Q) = mat(1/2, 0; 0, -1), quad
    vb(P)(t) = vb(Theta)(t) = mat(cos(t), sin(t); -sin(t), cos(t))
$
In either case $vb(R)$ is not unique: other branches of the logarithm give $vb(R) = op("diag")(1/2 + (2k + 1) j, -1 + (2m + 1) j)$ for $T = pi$ and $vb(R) = op("diag")(1/2 + k j, -1 + m j)$ for $T = 2 pi$ ($k, m$ integers), each with its own $vb(P)(t)$. For a given $T$, the Floquet multipliers ($op("eig")(vb(Q))$) are unique. For every choice, the real parts of the Floquet exponents (the eigenvalues of $vb(R)$) are $1/2$ and $-1$.

#emph[Stability:] In the coordinates $vbu(z) = vb(P)^(-1)(t) vbu(x)$ the dynamics are $dot(vbu(z)) = vb(R) vbu(z)$. A Lyapunov transformation preserves convergence and divergence. $vb(R)$ has an eigenvalue with $Re(r_1) = 1/2 > 0$, equivalently a Floquet multiplier with $abs(-e^(pi\/2)) approx 4.81 > 1$, so the system is unstable. This is consistent with HW02, where almost every solution grows like $e^(t\/2)$.
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

The system matrix is a scalar function times a fixed skew-symmetric matrix:
$
    vb(A)(t) = a(t) vb(S), quad a(t) = t sin(t) cos(sqrt(t)), quad vb(S) = mat(0, 1; -1, 0), quad vb(S)^TT = -vb(S)
$
Here $a(t)$ is real for every $t$: for $t < 0$, $cos(sqrt(t)) = cosh(sqrt(abs(t)))$. So the argument below holds for any $t_0$. Take $V(vbu(x)) = norm(vbu(x))_2^2$, i.e. $vb(P) = vb(I)$, so that $vb(Q)(t) = -[vb(A)^TT (t) + vb(A)(t)]$. Then
$
    dv(, t) norm(vbu(x)(t))_2^2
    = vbu(x)^TT (t) [vb(A)^TT (t) + vb(A)(t)] vbu(x)(t)
    = a(t) vbu(x)^TT (t) [vb(S)^TT + vb(S)] vbu(x)(t)
    = 0
$
so the norm never changes from its initial value:
$
    norm(vbu(x)(t))_2 = norm(vbu(x)_0)_2 = sqrt((-1)^2 + 2^2) = sqrt(5) approx 2.236 quad forall t >= t_0 quad qed
$

The solution itself is a rotation. Since $vb(A)(t) vb(A)(tau) = a(t) a(tau) vb(S)^2 = vb(A)(tau) vb(A)(t)$, the STM is the exponential of the integrated system matrix (the commuting case from #link("https://github.com/vaisriv/enae641-hw02")[HW02]). With $vb(S)^2 = -vb(I)$,
$
    vb(Phi)(t, t_0) = e^(theta vb(S)) = cos(theta) vb(I) + sin(theta) vb(S) = mat(cos(theta), sin(theta); -sin(theta), cos(theta)), quad
    theta(t, t_0) = integral_(t_0)^t sigma sin(sigma) cos(sqrt(sigma)) dd(sigma)
$
So $vbu(x)(t) = vb(Phi)(t, t_0) vbu(x)_0$ is $vbu(x)_0$ rotated by the angle $theta(t, t_0)$. This angle has no elementary closed form, but it is not needed for the norm.

#emph[Stability:] The argument above holds for every $vbu(x)_0$ and every $t_0$, so $norm(vbu(x)(t))_2 = norm(vbu(x)_0)_2$, i.e. $norm(vb(Phi)(t, t_0))_(i,2) = 1$:
- #emph[Stable isL:] $norm(vbu(x)(t))_2 <= gamma norm(vbu(x)_0)_2$ with $gamma = 1$, independent of $t_0$, so the stability is uniform. This matches the lecture test: $vb(P) = vb(I)$ gives $vb(Q) = vb(0)$, which is positive semi-definite.
- #emph[Not asymptotically stable:] $norm(vbu(x)(t))_2 = norm(vbu(x)_0)_2 arrow.not 0$ for any $vbu(x)_0 != vb(0)$, so it is not exponentially stable either.
- #emph[Not unstable:] every solution stays bounded.

The state keeps circling at radius $sqrt(5)$ about the origin. It never escapes and never decays:
$
    norm(vbu(x)(t))_2 = sqrt(5) quad => quad "stable isL (uniformly), but not asymptotically stable" quad qed
$
<hwk:s03>

#pagebreak(weak: true)

#prob[
    + For the F16 model in #link(<hwk:p01>)[Problem 1], find a constant, positive-definite matrix $vb(P)$ so that
        $
            dv((vbu(x)^TT (t) vb(P) vbu(x)(t)), t, style: "large") = -norm(vbu(x)(t))^2
        $
        (Use `MATLAB`'s ```matlab lyap()``` function, but read its documentation carefully!) Verify this $vb(P)$ is positive definite.
    <hwk:p04a>

    + Use this $(vb(P), vb(Q))$ pair to bound the exponential convergence of $||vbu(x)(t)||$ from an arbitrary $vbu(x)_0$.
    <hwk:p04b>
] <hwk:p04>

+
    For a constant $vb(P)$ and $dot(vbu(x)) = vb(A) vbu(x)$,
    $
        dv(, t) (vbu(x)^TT vb(P) vbu(x)) = dot(vbu(x))^TT vb(P) vbu(x) + vbu(x)^TT vb(P) dot(vbu(x)) = vbu(x)^TT (vb(A)^TT vb(P) + vb(P) vb(A)) vbu(x)
    $
    This equals $-norm(vbu(x))_2^2 = -vbu(x)^TT vb(I) vbu(x)$ for every $vbu(x)$ exactly when $vb(P)$ solves the algebraic Lyapunov equation with $vb(Q) = vb(I)$:
    $
        vb(A)^TT vb(P) + vb(P) vb(A) = -vb(I)
    $
    `MATLAB`'s ```matlab lyap(A, Q)``` solves $vb(A) vb(X) + vb(X) vb(A)^TT + vb(Q) = vb(0)$, which has the transpose on the other factor. Calling ```matlab lyap(A, eye(4))``` would therefore solve the wrong equation. That mistake would not show up in a positive-definiteness check either. Its solution is also positive definite, since $vb(A)^TT$ is Hurwitz too. Passing $vb(A)^TT$ instead, ```matlab P = lyap(A', eye(4))``` solves $vb(A)^TT vb(P) + vb(P) vb(A) + vb(I) = vb(0)$ and gives
    $
        vb(P) = mat(
            27.17, -3.394, 8.074, 0.2346;
            -3.394, 5.573, -16.59, -1.434;
            8.074, -16.59, 57.44, 4.671;
            0.2346, -1.434, 4.671, 0.5930;
        ) quad qed
    $
    #emph[Positive definiteness:] $vb(P)$ is symmetric, and all four of its eigenvalues are positive (a Cholesky factorization also succeeds). By Rayleigh--Ritz, $vbu(x)^TT vb(P) vbu(x) >= lambda_min (vb(P)) norm(vbu(x))_2^2 > 0$ for every $vbu(x) != vb(0)$:
    $
        lambda(vb(P)) = {0.1873, 0.6929, 25.20, 64.69} quad => quad lambda_min (vb(P)) = 0.1873 > 0 quad => quad vb(P) "is positive definite" quad qed
    $
    #raw(read("../outputs/text/s04a.txt"), block: true)
<hwk:s04a>

+
    Use the pair $(vb(P), vb(Q)) = (vb(P), vb(I))$ from #link(<hwk:s04a>)[Part A], with
    $
        rho = lambda_max (vb(P)) = 64.69, quad eta = lambda_min (vb(P)) = 0.1873, quad nu = lambda_min (vb(Q)) = 1
    $
    With $V = vbu(x)^TT vb(P) vbu(x)$, Rayleigh--Ritz gives $eta norm(vbu(x))_2^2 <= V <= rho norm(vbu(x))_2^2$. Since $nu > 0$,
    $
        dot(V) = -vbu(x)^TT vb(Q) vbu(x) <= -nu norm(vbu(x))_2^2 <= -nu/rho V
        quad => quad
        V(t) <= V(t_0) e^(-(nu\/rho)(t - t_0))
    $
    by the comparison lemma. Then $eta norm(vbu(x)(t))_2^2 <= V(t)$ and $V(t_0) <= rho norm(vbu(x)_0)_2^2$ give the lecture bound
    $
        norm(vbu(x)(t))_2 <= sqrt(rho/eta) e^(-lambda (t - t_0)) norm(vbu(x)_0)_2, quad
        gamma = sqrt(rho/eta) = 18.58, quad
        lambda = nu/(2 rho) = #qty(0.007729, "1/s")
    $
    which holds from an arbitrary $vbu(x)_0$ and $t_0$:
    $
        norm(vbu(x)(t))_2 <= 18.58 e^(-0.007729 (t - t_0)) norm(vbu(x)_0)_2 quad qed
    $
    This is a valid but conservative bound:
    - The guaranteed rate $lambda = #qty(0.00773, "1/s")$ is about half the true slowest decay rate $abs(sigma_1) = #qty(0.01425, "1/s")$ from #link(<hwk:s01a>)[Problem 1a]. Since $lambda = nu\/(2 rho)$, this gap depends on $rho$ alone, not on $eta$. $rho$ is set mainly by the phugoid. Restricted to the phugoid's plane, $vb(P)$ already has largest eigenvalue $approx 59.7$, above the $1\/(2 abs(sigma_1)) = 35.1$ that a mode with orthogonal, equal-length $Re{vbu(v)_1}$ and $Im{vbu(v)_1}$ would give. Here those vectors have lengths $1.03$ and $0.70$ and are #qty(77, "degree") apart, so the phugoid is not "round" in these coordinates.
    - The overshoot factor $gamma = 18.58$ is far above the true worst-case gain $max_t norm(e^(vb(A) t))_(i,2) = 1.127$, which is reached at $t approx #qty(17.5, "s")$. Even at the guaranteed rate, the smallest valid constant is $max_t norm(e^(vb(A) t))_(i,2) e^(lambda t) = 1.34$. $gamma = sqrt(rho\/eta)$ is large because the directions of $rho$ and $eta$ come from different modes. $rho$ comes mostly from the slow phugoid, while the eigenvector for $eta$ lies in the fast short-period plane, and their decay rates differ by a factor of about $abs(sigma_3)\/abs(sigma_1) approx 119$.
    - For comparison, let $vb(P)_m$ be the modal transformation of #link(<hwk:s01a>)[Problem 1a] (not the Lyapunov $vb(P)$ above). It gives the much tighter $norm(vbu(x)(t))_2 <= norm(vb(P)_m)_(i,2) norm(vb(P)_m^(-1))_(i,2) e^(sigma_1 t) norm(vbu(x)_0)_2 = 2.62 e^(-0.01425 t) norm(vbu(x)_0)_2$. The Lyapunov bound is looser, but it needs no eigenvectors, and the same argument carries over to a time-varying $vb(A)(t)$.
    #figure(
        image("../outputs/figures/s04b.png", width: 90%),
        caption: [Worst-case normalized response $norm(e^(vb(A) t))_(i,2) = max_(vbu(x)_0) norm(vbu(x)(t))_2 \/ norm(vbu(x)_0)_2$ against the Lyapunov bound $gamma e^(-lambda t)$.],
    )
    #raw(read("../outputs/text/s04b.txt"), block: true)
<hwk:s04b>

#pagebreak(weak: true)

#prob[
    + For the F16 model in #link(<hwk:p01>)[Problem 1], what are the eigenvalues of $vb(A) + vb(A)^TT$? Do these eigenvalues allow you to conclude the kind of exponential stability result that you obtained in #link(<hwk:p01>)[Problem 1] and #link(<hwk:p04>)[Problem 4]?
    <hwk:p05a>

    + Is there a contradiction between the results of #link(<hwk:p05a>)[Part A] and that of #link(<hwk:p01>)[Problem 1] and #link(<hwk:p04>)[Problem 4]? Explain carefully in light of Lyapunov theory.

        #smallcaps[Hint:] What $vb(P)$ are you implicitly using in the test of #link(<hwk:p05a>)[Part A]? Should you expect an arbitrary choice such as this to satisfy Lyapunov's requirement?
    <hwk:p05b>
] <hwk:p05>

+
    `MATLAB` gives
    $
        lambda(vb(A) + vb(A)^TT) = {-7.327, -0.2504, 0.04021, 0.7006} quad qed
    $
    These are the eigenvalues in the test from lecture: along solutions,
    $
        dv(, t) norm(vbu(x))_2^2 = vbu(x)^TT (vb(A)^TT + vb(A)) vbu(x)
    $
    The theorem says the system is stable isL if $lambda_max (vb(A)^TT + vb(A)) <= 0$. If instead $vb(A)^TT + vb(A) <= -epsilon vb(I)$, i.e. $vb(Q) = -(vb(A)^TT + vb(A))$ is uniformly positive definite, it is exponentially stable. Here two eigenvalues are positive, so $vb(Q)$ is indefinite. Neither condition holds, and the test cannot give exponential stability, or even stability isL. One entry already rules out the exponential test: $a_(33) = 0$, so $vbu(e)_3^TT (vb(A) + vb(A)^TT) vbu(e)_3 = 2 a_(33) = 0$, and $vb(A) + vb(A)^TT$ cannot be negative definite. Along the eigenvector $vbu(w)$ of the largest eigenvalue, $norm(vbu(x))_2^2$ actually starts out increasing: $dv(, t) norm(vbu(x))_2^2 = 0.7006 norm(vbu(x)_0)_2^2 > 0$ at $vbu(x)_0 = vbu(w)$ (#link(<hwk:s05b>)[Part B]).
    $
        lambda_max (vb(A) + vb(A)^TT) = 0.7006 > 0 quad => quad "no: the test is inconclusive" quad qed
    $
    #raw(read("../outputs/text/s05a.txt"), block: true)
<hwk:s05a>

+
    There is no contradiction. The test in #link(<hwk:s05a>)[Part A] is Lyapunov's test with the particular choice $vb(P) = vb(I)$: $V = vbu(x)^TT vb(I) vbu(x) = norm(vbu(x))_2^2$, and $vb(Q) = -(vb(A)^TT vb(I) + vb(I) vb(A)) = -(vb(A) + vb(A)^TT)$. Lyapunov's theorem is only a sufficient condition. When a particular $(vb(P), vb(Q))$ pair fails, the result is inconclusive, not a proof of instability, and a different $vb(P)$ may succeed.

    The converse theorem says what to expect. $vb(A)$ is Hurwitz (#link(<hwk:s01a>)[Problem 1a]), so for _any_ positive definite $vb(Q)$ there exists a positive definite $vb(P)$ with $vb(A)^TT vb(P) + vb(P) vb(A) = -vb(Q)$. #link(<hwk:s04a>)[Problem 4] found that $vb(P)$ for $vb(Q) = vb(I)$. The theorem does not say the reverse, that an arbitrary positive definite $vb(P)$ yields a positive definite $vb(Q)$. Choosing $vb(P) = vb(I)$ succeeds only if $vb(A) + vb(A)^TT$ is negative definite, i.e. the Euclidean norm $norm(vbu(x))_2$ shrinks along _every_ solution at _every_ instant. That is much stronger than all eigenvalues lying in the left half-plane, and this $vb(A)$ does not satisfy it. $vb(A)$ is far from normal: its couplings are lopsided (e.g. $a_(42) = -4.78$ but $a_(24) = 0.907$) and its states have very different scales. So the spherical level sets of $norm(vbu(x))_2^2$ are not adapted to the dynamics. The ellipsoidal level sets of $vbu(x)^TT vb(P) vbu(x)$ from #link(<hwk:s04a>)[Problem 4] are.

    Simulation shows that $norm(vbu(x))_2$ really does grow for a while from some initial states:
    - From $vbu(x)_0 = vbu(w)$ (the top eigenvector of $vb(A) + vb(A)^TT$), $norm(vbu(x))_2$ rises to $1.035 norm(vbu(x)_0)_2$ at $t approx #qty(0.16, "s")$.
    - From the worst-case unit initial state $vbu(x)_w$ for $t = #qty(17.5, "s")$, it rises to $1.127 norm(vbu(x)_0)_2$.

    So $norm(vbu(x))_2^2$ is not nonincreasing along all solutions, and no argument based on $vb(P) = vb(I)$ alone could certify the decay. Along the same solutions $V = vbu(x)^TT vb(P) vbu(x)$ decreases monotonically ($dot(V) = -norm(vbu(x))_2^2 < 0$). Both solutions are still decaying at the end of the plot, and they converge to zero, as #link(<hwk:s01a>)[Problem 1a] and #link(<hwk:s04b>)[Problem 4b] guarantee. An arbitrary $vb(P)$ such as $vb(I)$ should not be expected to satisfy Lyapunov's requirement. The part (a) test simply used the wrong $vb(P)$ for this $vb(A)$:
    $
        vb(P) = vb(I) => vb(Q) "indefinite (inconclusive)", quad vb(P) = "lyap"(vb(A)^TT, vb(I)) => vb(Q) = vb(I) succ 0 quad => quad "no contradiction" quad qed
    $
    #figure(
        image("../outputs/figures/s05b.png", width: 100%),
        caption: [Euclidean norm over #qty(60, "s") (left) and over the first second (middle), and $V = vbu(x)^TT vb(P) vbu(x)$ (right), each normalized by its initial value, from $vbu(x)_0 = vbu(w)$ and from the worst-case initial state $vbu(x)_w$.],
    )
    #raw(read("../outputs/text/s05b.txt"), block: true)
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

+
    The same bound is used in all four parts. Let $V(vbu(x), t) = vbu(x)^TT vb(P)(t) vbu(x)$, with $rho = max_t lambda_max (vb(P)(t)) < oo$ and $eta = min_t lambda_min (vb(P)(t)) > 0$. By Rayleigh--Ritz,
    $
        eta norm(vbu(x))_2^2 <= V <= rho norm(vbu(x))_2^2
    $
    and along solutions, by the Lyapunov equation $-vb(Q) = dot(vb(P)) + vb(A)^TT vb(P) + vb(P) vb(A)$,
    $
        dot(V) = -vbu(x)^TT vb(Q)(t) vbu(x) <= -nu(t) norm(vbu(x))_2^2
    $
    To compare $dot(V)$ with $V$, the sign of $nu(t)$ decides which side of the Rayleigh--Ritz bound applies:
    - Where $nu(t) >= 0$, use $norm(vbu(x))_2^2 >= V\/rho$: $quad -nu norm(vbu(x))_2^2 <= -(nu\/rho) V$.
    - Where $nu(t) < 0$, use $norm(vbu(x))_2^2 <= V\/eta$: $quad -nu norm(vbu(x))_2^2 = abs(nu) norm(vbu(x))_2^2 <= (abs(nu)\/eta) V$.
    Both cases are covered by
    $
        dot(V) <= -beta(t) V, quad
        beta(t) = min lr({(nu(t))/rho, (nu(t))/eta}) = (nu_(+) (t))/rho - (nu_(-) (t))/eta, quad
        nu_(plus.minus) = max lr({plus.minus nu, 0})
    $
    Comparing with $dot(w) = -beta(t) w$, $w(t_0) = V(t_0)$ (comparison lemma) gives $V(t) <= V(t_0) e^(-integral_(t_0)^t beta(sigma) dd(sigma))$. Then
    $
        norm(vbu(x)(t))_2 <= sqrt(rho/eta) exp(-1/2 integral_(t_0)^t beta(sigma) dd(sigma)) norm(vbu(x)_0)_2 quad (*)
    $
    so each case comes down to the behavior of $integral_(t_0)^t beta$.

    To show a conclusion cannot be strengthened, I use the scalar test system $dot(x) = -1/2 nu(t) x$ with $P = 1$ ($rho = eta = 1$). Its Lyapunov equation gives $Q = -2(-1/2 nu) = nu(t)$ exactly. Its solution $x(t) = x_0 exp(-1/2 integral_(t_0)^t nu)$ makes $(*)$ an equality, and it satisfies every hypothesis whenever $nu(t)$ is bounded. Since $vb(Q)$ must be bounded, take $t_0 > 0$ in parts (a) and (c) and $t_0 >= 0$ in part (b).

    Here $nu(t) = 1\/t > 0$ for all $t >= t_0 > 0$, so $beta = 1\/(rho t)$ and
    $
        integral_(t_0)^t beta(sigma) dd(sigma) = 1/rho ln(t/t_0)
        quad => quad
        norm(vbu(x)(t))_2 <= sqrt(rho/eta) (t_0/t)^(1\/(2 rho)) norm(vbu(x)_0)_2
    $
    - #emph[Stable isL:] since $nu >= 0$, $dot(V) <= 0$ and $norm(vbu(x)(t))_2 <= sqrt(rho\/eta) norm(vbu(x)_0)_2$, with $gamma = sqrt(rho\/eta)$ independent of $t_0$. This is the lecture's positive semi-definite case.
    - #emph[Asymptotically stable:] $(t_0\/t)^(1\/(2 rho)) -> 0$ as $t -> oo$, so $norm(vbu(x)(t))_2 -> 0$ for every $vbu(x)_0$.
    - #emph[Uniformity cannot be concluded:] the bound drops below $epsilon$ only once $t >= t_0 (sqrt(rho\/eta)\/epsilon)^(2 rho)$, so the $T(epsilon, t_0)$ it certifies grows with $t_0$. Some systems that meet the hypotheses really are not UAS. The test system $dot(x) = -x\/(2 t)$ has $nu = 1\/t$ and $x(t) = x_0 sqrt(t_0\/t)$, so it needs $T(epsilon, t_0) = t_0 (epsilon^(-2) - 1)$, which grows with $t_0$, like the lecture example $vb(Phi)(t, t_0) = (1 + t_0)\/(1 + t)$.
    - #emph[Not exponentially:] $vb(Q)(t)$ is positive definite at every instant but not uniformly ($inf_t nu(t) = 0$), so the lecture ES test does not apply, and the decay above is only algebraic. The same test system decays like $x_0 sqrt(t_0\/t)$, which no envelope $gamma e^(-lambda (t - t_0))$ can bound as $t -> oo$.
    $
        nu(t) = 1/t quad => quad "asymptotically stable (UAS and ES cannot be concluded)" quad qed
    $
<hwk:s06a>

+
    Here $nu(t) = -c e^(-t) < 0$ for all $t$, so $vb(Q)(t)$ is never positive semi-definite and neither lecture test applies. Every instant falls in the $nu < 0$ case, so $beta = -(c\/eta) e^(-t)$ and
    $
        integral_(t_0)^t beta(sigma) dd(sigma) = -c/eta (e^(-t_0) - e^(-t)) >= -c/eta e^(-t_0)
    $
    The integral is bounded below, so by $(*)$, for all $t >= t_0 >= 0$,
    $
        norm(vbu(x)(t))_2 <= sqrt(rho/eta) exp(c/(2 eta) (e^(-t_0) - e^(-t))) norm(vbu(x)_0)_2 <= sqrt(rho/eta) e^(c\/(2 eta)) norm(vbu(x)_0)_2
    $
    This is stability isL, with $gamma = sqrt(rho\/eta) e^(c\/(2 eta))$ independent of $t_0$. $V$ may grow, but the total destabilizing contribution $integral_(t_0)^oo abs(nu(sigma)) dd(sigma) = c e^(-t_0)$ is finite, so it grows by at most a fixed factor.

    Asymptotic stability cannot be concluded. The test system $dot(x) = (c\/2) e^(-t) x$ has $nu(t) = -c e^(-t)$, yet
    $
        x(t) = x_0 exp(c/2 (e^(-t_0) - e^(-t))) -> x_0 e^(c e^(-t_0) \/ 2) != 0
    $
    $
        nu(t) = -c e^(-t) quad => quad "stable isL (uniformly)" quad qed
    $
<hwk:s06b>

+
    Here $nu(t) = 1\/t^2 > 0$, so $vb(Q)(t)$ is positive definite at every instant, $dot(V) <= 0$, and (as in the lecture theorem) $norm(vbu(x)(t))_2 <= sqrt(rho\/eta) norm(vbu(x)_0)_2$: stable isL, with $gamma$ independent of $t_0$. However, $nu$ decays fast enough that its integral converges:
    $
        integral_(t_0)^t beta(sigma) dd(sigma) = 1/rho (1/t_0 - 1/t) -> 1/(rho t_0) < oo
    $
    so $(*)$ only gives $norm(vbu(x)(t))_2 <= sqrt(rho\/eta) e^(-(1\/t_0 - 1\/t)\/(2 rho)) norm(vbu(x)_0)_2$, which does not tend to zero. This is a real limitation, not a weak bound. The test system $dot(x) = -x\/(2 t^2)$ has $nu = 1\/t^2$, and
    $
        x(t) = x_0 exp(-1/2 (1/t_0 - 1/t)) -> x_0 e^(-1\/(2 t_0)) != 0
    $
    So asymptotic stability cannot be concluded. A $vb(Q)(t)$ that is positive definite at every instant is not enough when the total dissipation $integral_(t_0)^oo nu$ is finite.
    $
        nu(t) = 1/t^2 quad => quad "stable isL (uniformly)" quad qed
    $
<hwk:s06c>

+
    Here $nu(t) = 2 - 3 sin(t)$ ranges over $[-1, 5]$. It is negative when $sin(t) > 2\/3$, i.e. on $(t_1, pi - t_1) + 2 k pi$ with $t_1 = arcsin(2\/3)$. So $vb(Q)(t)$ is not positive semi-definite for part of every period, and neither lecture test applies. Still, $nu$ is positive on average. Over one period, using $cos(t_1) = sqrt(5)\/3$:
    $
        I_(-) & = integral_0^(2 pi) nu_(-) dd(t) = integral_(t_1)^(pi - t_1) (3 sin(t) - 2) dd(t) = 2 sqrt(5) - 2 (pi - 2 t_1) approx 1.108 \
        I_(+) & = integral_0^(2 pi) nu_(+) dd(t) = integral_0^(2 pi) nu dd(t) + I_(-) = 4 pi + I_(-) approx 13.674
    $
    $beta$ is $2 pi$-periodic, so its integral over any window of length $2 pi$ is
    $
        mu = integral_t^(t + 2 pi) beta(sigma) dd(sigma) = I_(+)/rho - I_(-)/eta
    $
    Write $t - t_0 = 2 pi k + r$ with $0 <= r < 2 pi$. The $k$ full periods contribute $k mu$, and the leftover piece contributes at least $-I_(-)\/eta$. So if $mu > 0$,
    $
        integral_(t_0)^t beta(sigma) dd(sigma) >= mu (t - t_0)/(2 pi) - mu - I_(-)/eta
        quad => quad
        norm(vbu(x)(t))_2 <= sqrt(rho/eta) e^((mu + I_(-)\/eta)\/2) e^(-mu (t - t_0)\/(4 pi)) norm(vbu(x)_0)_2
    $
    Both constants are independent of $t_0$, so the system is uniformly exponentially stable with $lambda = mu\/(4 pi)$ whenever
    $
        mu > 0 quad <=> quad rho/eta < kappa^* = I_(+)/I_(-) = 1 + (4 pi)/(2 sqrt(5) - 2 pi + 4 arcsin(2\/3)) approx 12.343
    $
    In particular this always holds when $vb(P)$ is a constant multiple of $vb(I)$ ($rho = eta$, so $mu = 4 pi\/rho$). In the boundary case $rho\/eta = kappa^*$, $mu = 0$ and the same window argument gives $integral_(t_0)^t beta >= -I_(-)\/eta$. By $(*)$ that is uniform stability isL, with $gamma = sqrt(rho\/eta) e^(I_(-)\/(2 eta))$. The bound $(*)$ does not decay there, so AS cannot be concluded from this pair.

    If $rho\/eta > kappa^*$, then $mu < 0$, the bound grows, and the test is inconclusive. This is not just a loose bound. Take $vb(P) = op("diag")(rho, eta)$, $vb(Q)(t) = nu(t) vb(I)$, and $vb(A)(t) = vb(P)^(-1) [-1/2 nu(t) vb(I) + k(t) vb(S)]$ with $vb(S) = mat(0, 1; -1, 0)$. For any bounded scalar $k(t)$, these satisfy the Lyapunov equation with $min_k lambda_k (vb(Q)(t)) = nu(t)$. In $vbu(y) = vb(P)^(1\/2) vbu(x)$ the dynamics are $dot(vbu(y)) = [-1/2 nu vb(P)^(-1) + k\/sqrt(rho eta) thin vb(S)] vbu(y)$, so the $k(t)$ term only rotates $vbu(y)$ and leaves $V = norm(vbu(y))_2^2$ unchanged. Start at $t_0 = 0$ with $vbu(x)_0$ on the $rho$-axis ($nu(0) = 2 > 0$). Let $k(t) = plus.minus (pi\/2) sqrt(rho eta)\/delta$ on windows of length $delta$ at each zero crossing of $nu$, and $0$ elsewhere, which makes a quarter turn at each crossing. Then the state lies along the $eta$-axis while $nu < 0$, where $V$ grows at rate $abs(nu)\/eta$. It lies along the $rho$-axis while $nu > 0$, where $V$ decays only at rate $nu\/rho$. Each period multiplies $V$ by $e^(-mu)$, up to an error that vanishes as $delta -> 0$. For any $rho\/eta > kappa^*$ a small enough $delta$ makes that factor exceed $1$, so this system is unstable.
    $
        nu(t) = 2 - 3 sin(t) quad => quad cases(
            "uniformly exponentially stable" & "if" rho\/eta < kappa^* approx 12.343,
            "stable isL (uniformly)" & "if" rho\/eta = kappa^*,
            "inconclusive" & "if" rho\/eta > kappa^*,
        ) quad qed
    $
<hwk:s06d>

#pagebreak(weak: true)

#prob[
    Generalize your observations from #link(<hwk:p06>)[Problem 6]. What general conditions on $nu(t)$ would allow you to conclude isL stability? Asymptotic stability? Exponential stability?

    #smallcaps[Note:] This question is asking you to go beyond the comparatively straightforward results discussed in lecture; for example to demonstrate formally that stability can sometimes be concluded if $nu(t) < 0$ for some or even all $t >= t_0$, etc.
] <hwk:p07>

All four cases of #link(<hwk:p06>)[Problem 6] follow from one bound. Assume $vb(P)(t)$ is bounded and uniformly positive definite ($eta vb(I) <= vb(P)(t) <= rho vb(I)$) and $vb(Q)(t)$ is bounded ($norm(vb(Q)(t))_i <= q$, so $nu(t) >= -q$). With $V = vbu(x)^TT vb(P) vbu(x)$, the argument at the start of #link(<hwk:s06a>)[Problem 6a] gives
$
    dot(V) <= -beta(t) V, quad beta(t) = (nu_(+) (t))/rho - (nu_(-) (t))/eta
    quad => quad
    norm(vbu(x)(t))_2 <= sqrt(rho/eta) exp(-1/2 B(t, t_0)) norm(vbu(x)_0)_2, quad
    B(t, t_0) = integral_(t_0)^t beta(sigma) dd(sigma)
$
where $nu_(plus.minus) = max{plus.minus nu, 0}$. The positive part of $nu$ is weighted by $1\/rho$ and the negative part by $1\/eta$. Where $nu < 0$, the state can sit in the direction where $vb(P)$ is smallest, making $norm(vbu(x))_2^2$ as large as $V\/eta$. So what this $(vb(P), vb(Q))$ pair can certify is decided by the running integral $B(t, t_0)$, not by the sign of $nu$ at each instant:

- #emph[Stable isL] if $B$ is bounded below: $B(t, t_0) >= -M(t_0) > -oo$ for all $t >= t_0$. Then $gamma(t_0) = sqrt(rho\/eta) e^(M(t_0)\/2)$, and the stability is uniform if $M$ does not depend on $t_0$. A sufficient condition on $nu$ alone is integrability of its negative part:
    $
        integral_(t_0)^oo nu_(-) (sigma) dd(sigma) < oo
        quad => quad
        norm(vbu(x)(t))_2 <= sqrt(rho/eta) exp(1/(2 eta) integral_(t_0)^oo nu_(-) (sigma) dd(sigma)) norm(vbu(x)_0)_2
    $
    since $B(t, t_0) >= -1/eta integral_(t_0)^t nu_(-)$. This stability is automatically uniform over all $t_0 >= t_s$ (for any fixed earliest start time $t_s$), because $integral_(t_0)^oo nu_(-) <= integral_(t_s)^oo nu_(-)$. So stability can be concluded even if $nu(t) < 0$ for some or all $t >= t_0$, as long as the total negative contribution is finite (#link(<hwk:s06b>)[6b]). The lecture condition, $vb(Q)$ positive semi-definite ($nu >= 0$), is the special case $nu_(-) equiv 0$, $M = 0$. If $nu(t) < 0$ for all $t$, then $B <= 0$, so stability isL is the most this pair can certify; AS and ES cannot be concluded.

- #emph[Asymptotically stable] if $B(t, t_0) -> +oo$ as $t -> oo$, since then $norm(vbu(x)(t))_2 -> 0$ for every $vbu(x)_0$. This already includes stability isL: $B$ is continuous with $B(t_0, t_0) = 0$, so a $B$ that tends to $+oo$ is bounded below. A sufficient condition on $nu$ alone is
    $
        integral_(t_0)^oo nu_(-) (sigma) dd(sigma) < oo quad "and" quad integral_(t_0)^oo nu_(+) (sigma) dd(sigma) = oo
    $
    The divergence matters. $nu = 1\/t$ (#link(<hwk:s06a>)[6a]) and $nu = 1\/t^2$ (#link(<hwk:s06c>)[6c]) are both positive, but only the first has a divergent integral, and only the first gives AS. If $B$ diverges more slowly for later $t_0$ (as $ln(t\/t_0)$ does), the time the bound needs to reach a given $epsilon$ grows with $t_0$, so uniformity cannot be concluded from this bound.

- #emph[Exponentially stable] if $B$ grows at least linearly,
    $
        B(t, t_0) >= 2 lambda (t - t_0) - M quad forall t >= t_0, quad "for some" lambda > 0, M < oo
    $
    which gives $norm(vbu(x)(t))_2 <= sqrt(rho\/eta) e^(M\/2) e^(-lambda (t - t_0)) norm(vbu(x)_0)_2$. If $lambda$ and $M$ do not depend on $t_0$, it is UES. A convenient sufficient condition is a uniformly positive windowed integral: for some $T > 0$ and $mu > 0$,
    $
        integral_t^(t + T) beta(sigma) dd(sigma) = integral_t^(t + T) [(nu_(+) (sigma))/rho - (nu_(-) (sigma))/eta] dd(sigma) >= mu quad forall t >= t_0
    $
    Chop $[t_0, t]$ into windows of length $T$, as in the lecture proof that UAS implies UES. Each full window contributes at least $mu$, and the leftover piece at least $-q T\/eta$, so $lambda = mu\/(2 T)$ and $M = mu + q T\/eta$, independent of $t_0$: UES. The lecture case ($vb(Q)$ uniformly positive definite, $nu(t) >= nu_min > 0$) is the special case $mu = nu_min T\/rho$, which recovers $lambda = nu_min\/(2 rho)$. #link(<hwk:s06d>)[6d] satisfies the condition with $nu < 0$ on part of every period, provided $rho\/eta < kappa^* approx 12.343$. A version that uses $nu$ alone and works for every $rho\/eta$: if $integral_(t_0)^oo nu_(-) = N < oo$ and $integral_t^(t + T) nu_(+) >= mu' > 0$ for all $t$, then $B(t, t_0) >= (mu'\/rho)((t - t_0)\/T - 1) - N\/eta$. Since $N <= integral_(t_s)^oo nu_(-)$ for every $t_0 >= t_s$, the constant $M = mu'\/rho + N\/eta$ is uniform, and this is again UES.

These conditions are sufficient, not necessary. When they fail, the result is inconclusive. The bound comes from one particular $(vb(P), vb(Q))$ pair, and by the converse theorem a UES system has infinitely many other pairs, some of which may pass. The conditions also cannot be weakened using only $nu(t)$. For the scalar systems $dot(x) = a(t) x$ with constant $P = p$ (so $rho = eta = p$), $dot(V) = -nu x^2 = -(nu\/p) V$ holds with equality. The bound is then exact, and each condition on $B$ is also necessary for that family, as the test systems in #link(<hwk:p06>)[Problem 6] showed. The different weights $1\/rho$ and $1\/eta$ matter only when $nu$ changes sign, which is why the conditioning $rho\/eta$ of $vb(P)$ entered #link(<hwk:s06d>)[6d]. The rotating construction there shows that the $1\/eta$ weight on $nu_(-)$ cannot be improved when $rho != eta$. In summary, with $B(t, t_0) = integral_(t_0)^t [nu_(+)\/rho - nu_(-)\/eta] dd(sigma)$:
$
    "stable isL" & <== quad inf_(t >= t_0) B(t, t_0) > -oo quad (<== integral_(t_0)^oo nu_(-) < oo) \
    "asymptotically stable" & <== quad B(t, t_0) -> oo \
    "exponentially stable" & <== quad B(t, t_0) >= 2 lambda (t - t_0) - M, quad "uniformly (UES) if" integral_t^(t + T) beta >= mu > 0 thick forall t quad qed
$
<hwk:s07>

#pagebreak(weak: true)

== Code

#codly(header: [./src/index.m])
#raw(read("../src/index.m"), block: true, lang: "matlab") <code:index.m>
