# CoilRailSim Electromagnetic Acceleration Simulator

The following paper presents a physical derivation for the CoilRailSim application, as well as an algorithm for direct implementation and serves for educational purposes only. 
This algorithms simulates the influence of a solenoid onto a target in it's vicinity with all necessary physical parameters. 


## 1. Background

## 1.1. Magnetism of a Loop

To use coils for magnetic acceleration we need to determine the time-resolved fields and resulting forces which act upon a proximate probe in order to describe it's motion. For simplicity we approximate the coil as infinitely thin. A physical coil has $N$ turns of wire, length $l$ and usually a non-zero thickness $l\rightarrow0$. The reason for this idealization is the Biot-Savart law which maps the spatial propagation of magnetic flux in the proximity of flowing charge $I$ along a contour (or wire) $C$. More precisely, the law allows to calculate the resultant magnetic flux density $\vec{B}$ in 3-dimensional distance $\vec{x}$ from the loops center $\vec{x}_0$, given that the contour $C$ and the current $I$ are known.

$$ \vec{B}(\vec{x}) = \frac{\mu_0}{4\pi} \int_C \frac{I d\vec{s} \times \hat{r}'}{ \left | \vec{r}' \right |^2}$$

Here $\hat{r}'$ represents the unit vector in direction of $\vec{r}'$, $d\vec{s}$ is the infinitesimal segment along the contour $C$ in which the charges flow.

Coils are usually chosen in a concentric shape which automatically determines the contour $C$ (a circle) and allows analytical integration whose result can be written in closed-form. Hence, the solution for the resulting flux density of a short coil (crudely approximated as a loop) along the central perpendicular axis (in this case called "x" axis) yields the analytic solution

$$\vec{B}_\text{coil}(x) = \frac{\mu_0 NI R^2}{2(x^2+R^2)^{3/2}} \hat{x}$$

where $x$ is the scalar distance along the perpendicular axis. This formula is sub-optimal for our task, as we will explore in the upcoming section, but suffices for the current argument.

Ferromagnetic responses happen due to magnetization under a given external field $H$ which in our case will be the coil. The magnetic field which spans from the center of the loop along the $x$-axis is given by


$$\vec{H}(x) = \frac{\vec{B}_\text{coil}(x)}{\mu_0} = \frac{NI R^2}{2(x^2+R^2)^{3/2}} \hat{x}$$

where we assume that the vacuum permeability $\mu_0$ is scalar due to isotropy of the medium. $\vec{H}$ permeates every body in proximity and induces a magnetization in them based on their response function. To have the appropriate relation we have to determine the susceptibility $\chi_m$ which analogous to the permeability that relates $B$ and $H$, for linear materials (such as iron) acts as the linear relator between the $H$ and $M$-field as

$$ \vec{M} = \chi_m \vec{H}.$$

This is a linear reponse function for which we can compute the specific value for iron as

$$ \chi_\text{iron} = \mu_r-1 = \frac{\mu_\text{iron}}{\mu_0} - 1 \approx 5000.$$

The $M$-field is defined as the local magnetization gradient on an extended (magnetizable) manifold which is generalized for non-homogenous bodies as $M = \frac{d\mathcal{M}}{dV}$ where $\mathcal{M}$ is the total integrated magnetic moment. Since we can assume that the coin is homogeneous (at least in its ferromagnetic core) the integral simplifies to a linear extensive relation


$$ \vec{\mathcal{M}} = V \cdot \vec{M}   $$

i. e. the larger the body (either by volume $V$ or equivalently by mass $m=\rho V$) the stronger is it's magnetic dipole and response to external fields $H$ which relates to all forces. For the iron body the total magnetic dipole which builds up in vicinity to the coil would be

$$ \vec{\mathcal{M}}(x) = \chi_\text{iron}  V   \vec{H}(x) $$

which already shows that the magnetization within the iron body changes with distance to the coil - it increases with shrinking distance $x$. 

Every loop becomes an inductor as the field lines $\vec{H}(x)$ create a flux which demands magnetic energy or work. To have some relation between the current which passes through the loop and the resulting flux $\Phi$ we infer the "inductance" $\mathcal{L}$ which relates them linearly

$$\Phi = \mathcal{L} \cdot I$$

simultaneously it governs Faraday's induction law as it relates the (self-)induced voltage to the change in current

$$U_i = - \frac{d\Phi}{dt} = -\mathcal{L} \cdot \frac{dI}{dt}$$

where the negative sign represents Lenz's law. This is the ideal case definition where the inductance is constant as it can be obtained from the solenoid geometry and permeability

$$\mathcal{L} = \mu_{r} \mu_0 \frac{N^2 A}{l}$$

where $\mu_0$ is the vacuum permeability and $\mu_r$ the relative permeability of the magnetic core. Note that $\mu_r=1$ for a coil in vacuum. 

## 1.2. Solenoid Correction

In the previous part we found the magnetic potential by approximating a "thin" and "short" coil to a wire loop. 

$$\vec{H}_\text{loop}(x) = \frac{NI R^2}{2(x^2+R^2)^{3/2}} \hat{x}$$

Instead of absorbing all $N$ turns into an infinitely short loop which stands at $x=x_0 = 0$, we can reach better precision by using an integral form. In principle we integrate multiple loops with $N$ neighboring turns, each representing a separate loop, along the x-axis whose magnetic potentials interfere constructively. Since we have the 1-dimensional solution for a loop, for a coil we can integrate the loops discretely along the whole length $l$, such that the loop potentials at an arbitrary point $x$ accumulate. Each new loop will generate a shifted field by $l/N$ to it's neighbor as

$$ \vec{H}_\text{sol} (x) = \sum_{n=1}^{N} H_\text{loop}(x-x_0+\frac{n}{N}l) = \sum_{n=1}^{N} \frac{IR^2}{2( R^2 + (x - x_0 + l \cdot n/N)^2 )^{3/2}}$$

where $x$ is the arbitrary position of the projectile, $x_0$ is the start of the coil and $x_0 + \frac{n}{N}l$ the position of the $n$-th turn, until the end of the coil at $x_0+l$ is reached. 

Despite it's accuracy, the formula needs to be re-evaluated for every measurement of $H$, instead we can transfer the problem into continuous space to find a closed form

$$\vec{H}_\text{sol} (x) = \int_{x_0}^{x_0+l} \frac{NIR^2}{2l( R^2 + (x - dx')^2)} dx' = \frac{NI}{2l} \left [ \frac{x_0 + l - x}{\sqrt{R^2+(x_0+l-x)^2}} - \frac{x_0 - x}{\sqrt{R^2+(x_0-x)^2}} \right ]$$

which requires a single computation operation, instead of $N$. We obtain the final form, which is the standard *finite coil* potential

$$\boxed{\vec{H}_\text{sol} (x; x_0,N,l,I,R) = \frac{NI}{2l} \left [ \frac{x_0 + l - x}{\sqrt{R^2+(x_0+l-x)^2}} - \frac{x_0 - x}{\sqrt{R^2+(x_0-x)^2}} \right ] }$$

## 2. Dynamics

## 2.1. Biot-Savart Ansatz

Every force in classical mechanics can be derived from the potential energy gradient. Forces arise from slopes (gradients) in the corresponding scalar potential $U(x)$

$$\vec{F} = -\vec{\nabla} U.$$

The potential energy of a magnetic moment $\vec{\mathcal{M}}$ in an external field $\vec{B}$ is given by $U=-\vec{\mathcal{M}}(x) \cdot \vec{B}(x)$. We obtain the total force acting on a magnetized iron body at distance $x$ 


$$ \vec{F}(x) =  \frac{\partial}{\partial x} \left [  \vec{\mathcal{M}}(x) \cdot \vec{B}(x)  \right]$$

Since we know that $\vec{\mathcal{M}}$ is directly proportional to $\vec{B}$ at all times as $\vec{\mathcal{M}} = \frac{\chi V}{\mu} \vec{B}(x)$ we get

$$ 
\Rightarrow \vec{F}(x) = \frac{\chi  V}{\mu}  \frac{\partial}{\partial x} \left [  \vec{B}^2(x)  \right] = \chi \mu  V  \frac{\partial}{\partial x} \left [  \vec{H}_{sol}^2(x)  \right] 
$$

where $\mu$ is the net permeability of the medium $\mu = \mu_0 \mu_r$ - usually for air, $V$ is the projectile volume, and $\chi$ is the rel. susceptibility of the projectile material. Inserting the solution for $H_{sol}$ yields


$$ 
\Rightarrow \vec{F}(x) =  \frac{\chi \mu  V N^2 I^2 }{2l^2} 
\left ( 
\frac{(l-x+x_0)^2}{ ((l-x+x_0)^2+R^2)^{3/2} }
-\frac{1}{ \sqrt{(l-x+x_0)^2+R^2} }
-\frac{(x_0-x)^2}{ ((x-x_0)^2+R^2)^{3/2} }
+\frac{1}{ \sqrt{(x - x_0)^2+R^2} }
\right)
\left (
\frac{l-x+x_0}{ \sqrt{(l-x+x_0)^2+R^2}}
-\frac{x-x_0}{ \sqrt{(x-x_0)^2+R^2} } 
\right)
\hat{x}
$$

Since the force term is proportional to  $\chi \mu$ and large susceptibilities $\chi >> 1$ can be simplified to $\chi = \mu_r - 1 \approx \mu_r$ and if the projectile enters the solenoid it increases the solenoid permeability where $F \propto \mu_r \mu \rightarrow \mu_r^2 \mu_0$. This can lead to a blow-up if we assume uncoupled behaviour between projectile and solenoid as for iron $\mu_r \approx 5000$ leading to a false assumption for the force amplification. It cannot be assumed that the magnetization induced by the external field would lead to an identical magnet in the projectile. Instead, the projectile magnetization leads to a dipole moment which counter-acts the external solenoid field in it's vicinity which leads to a massive drop in the susceptibility. To correct for this we need to reduce the susceptibility to an "effective" value 

$$\chi \rightarrow\chi_\text{eff} = \frac{\chi}{1+\gamma \chi} \approx \frac{\chi}{\gamma \chi} = \frac{1}{\gamma}$$

When $\chi >> 1$, the physical property of the material ($\chi$) cancels out completely, leaving the effective susceptibility to be determined entirely by the probe geometry $\gamma$. The [demagnetizing factor](https://en.wikipedia.org/wiki/Demagnetizing_field#Demagnetizing_factor) $\gamma$ is a geometric tensor component that always sums to 1 across the three spatial dimensions ($\gamma_x + \gamma_y + \gamma_z = 1$). Assuming the external field is aligned along the long axis (x-axis) of the projectile, the table below shows the corresponding $\gamma$ values for different shapes:

|Shape|$\gamma$|
|-|-|
|Spherical|$\frac{1}{3}$|
|Cylinder of length $\mathcal{l}$ and radius $r$ (parallel)|$\frac{1}{1 + \frac{2\mathcal{l}}{r \sqrt{\pi}}}$|
|Disk of width $\mathcal{l}$ and radius $r$ (parallel)|$\frac{\pi \mathcal{l}}{8r}$|

Once the demagnetizing factor was determined, it can be plugged into the force formula and also exchange the projectile volume for the more accessible mass and density via $V = m / \rho$. The coil is placed to start at origin i. e. $x_0=0$ (stretching into positive direction) then we obtain the final force term which acts from the coil on the projectile

$$ 
\boxed {
\Rightarrow \vec{F}(x; \gamma, \mu, m, \rho, N,l,I,R) = \frac{\mu  m N^2 I^2 }{2\gamma\rho l^2}
\left ( 
\frac{(l-x)^2}{ ((l-x)^2+R^2)^{3/2} }
-\frac{1}{ \sqrt{(l-x)^2+R^2} }
-\frac{x^2}{ (x^2+R^2)^{3/2} }
+\frac{1}{ \sqrt{x^2+R^2} }
\right) 
\left ( 
\frac{l-x}{ \sqrt{(l-x)^2+R^2}}
+\frac{x}{ \sqrt{x^2+R^2} }
\right) \hat{x}
} 
$$


|Symbols|Quantity|
|-|-|
|$x$|Arbitrary position to coil along it's central axis|
|$\gamma$|Demagnetizing factor|
|$\chi$|Susceptibility of the projectile|
|$\mu$|Permeability of the medium|
|$m$|Mass of the projectile|
|$\rho$|Density of the projectile|
|$N$|Number of wire turns in the coil|
|$l$|Coil length|
|$I$|Current passing through the coil|
|$R$|Radius of the coil.|


## 2.2. Inductance Gradient Ansatz

Another ansatz which avoids heavy geometric calculations and tracking of parameters exploits that $U$ also equals the magnetic potential stored in the coil. Changes in the magnetic energy manifest in mechanical work i.e. accelerating the projectile which drains the coil energy stronger the closer it gets. Revisiting the definition of the reluctant force leads to

$$F(x) = \nabla W_\text{mag} = \nabla \left ( \frac{1}{2} \mathcal{L}(x) \cdot I^2 \right ) = -\frac{1}{2} \mathcal{L}'(x) \cdot I^2$$

which is a short elegant description and directly applicable with the caveat that the response $\mathcal{L}'(x)$ (i.e. the shape of $\mathcal{L}$ in dep. to projectile position) must be known. In trade off for the computational ease we loose accuracy. Note that the negative sign is not necessary as the potential is not conservative. The force is sensitive to fitting the correct inductance shape, which is merely approximated. Thus this ansatz alone sacrifices the fine, localized field structures that a full Biot-Savart integration yields. Nonetheless, the approach is interesting as all the previous geometry and material constants such as susceptibility vanish and the projectile motion is being purely governed by the inductance profile - similar to a potential valley governing the equations of motion.




## 3. Perturbation

Instead of using the inductance gradient to calculate the main force contribution which might introduce biases, we can use it to model the perturbation contributions.
Let's revisit the magnetic flux formula from section 1.1.

$$\Phi = \frac{\mathcal{L}}{N}\cdot I.$$

The self-induced voltage building up in the coil is given by the change of magnetic flux in time and always counteracts it's cause (Lenz's law)

$$V_L = -N \frac{d \Phi}{ dt} = - \frac{d (\mathcal{L} \cdot I)}{dt} = - \left ( \frac{d\mathcal{L}}{dt}I+\mathcal{L}\frac{dI}{ dt} \right )$$

since $\mathcal{L}$ changes with the projectile position (canonical variable which depends on time), we need to expand the composite dependence using the chain rule

$$\Rightarrow - V_L = \mathcal{L}\frac{dI}{ dt} + \frac{d\mathcal{L}}{dx} \frac{dx}{dt} I = \mathcal{L} \dot{I} + \mathcal{L}'Iv$$

where the first term is just the inductance definition $U=\mathcal{L}\dot{I}$ from Faraday's law of induction. The second term comes from the fact that $L$ is not constant, it varies due to moving projectile perturbing the field as it gets closer to the coil center. Hence, the spatial derivative does not describe the shape of $L$ along the coil axis but instead the inductive response to the change of the projectile position i.e. if the projectile is close to the coil center its movement will impact the inductance stronger than when its far from the coil. When the projectile comes near the coil it starts to raise the coil inductivity $\mathcal{L}$ as it provides a more favorable path for the field lines. This will reduce the current and lead to a damping effect - *Motional EMF*. Since the second part scales directly with velocity $v$, a fast flying projectile can make this term become significant where perturbation becomes comparable to the expected ideal values. 

## 3.1. First-Order Inductance Perturation Model

This section investigates the designing and modelling of the inductance profile, but only for illustrative purposes. For the actual calculation we will use a more elegant and phys. precise alternative which will be discussed in the next section.

The perturbation demands an inductance gradient for which we design a model in dependence to $x$. We start with the simple initial case in which the coil is empty and the projectile is far away. We can state that $\mathcal{L}$ is the conductance in air as the permeable coil medium i.e. $\mathcal{L}_\text{air}$. In the second extreme we fill the complete coil volume with the projectile material (i.e. steel, iron etc.) which shall raise the permeability significantly. This results in a maximum achievable inductance $\mathcal{L}_\text{max}$.
During the whole experiment the measured inductance is ensured to operate within these boundaries $\mathcal{L}_\text{air}  \leq \mathcal{L} \leq \mathcal{L}_\text{max}$. The both boundaries can be expressed as

$$\mathcal{L}_\text{air} = \mu_0 \frac{N^2 A}{l}$$

and for the filled case we can exploit the relative susceptibility $\chi$ of the projectile material.

$$\mathcal{L}_\text{max} = (1+\chi) \cdot \mathcal{L}_\text{air} $$

Next we take into account that there is significant mag. field leakage as the projectile does not fully fill up the coil volume and thus catches the flux only partially. Thus the change in inductance has to be corrected by relating the sizes of projectile and coil. Within the solenoid the projectile area with a perpendicular component to the field lines matter more than its length. Hence we relate the cross-sectional areas of solenoid and projectile. The inductance difference is multiplied with an *areal flux occupation* factor $\lambda$ to correct for the leakage

$$\lambda = \text{min} \left\{1, \frac{A_\text{projectile}}{A_\text{coil}} \right\}  = \text{min} \left\{1, \frac{A_\text{p}}{\pi R^2 } \right\} = \text{const.} \le 1$$

By further assuming that the projectile is fully contained within the solenoid we can exploit the homogeneity of it's internal magnetic field $H$. As a long ferromagnetic projectile moves deep inside the core, it fills the coil's volume uniformly, meaning the total magnetic flux and inductance reach a flat, maximum plateau. This introduces translational symmetry as no matter where the projectile is located within the coil it should distort the field lines equally. The inductance only changes significantly at the edges (the entry and exit zones) where the projectile transitions into and out of the magnetic field. To simulate this behaviour we need to design a plateau function $\Pi(x, x_0, L, \sigma)$ which is mainly constant in the center with smooth edges (whose transition width is characterized by $\sigma$). A neat way to achieve this is the double-sided sigmoid

$$\Pi(x) = \frac{1}{1+\text{exp} \left ( - \frac{x}{\sigma} \right ) } - \frac{1}{1+\text{exp} \left (  -\frac{x-l}{\sigma} \right ) } \equiv S(x/\sigma) - S((x-l) / \sigma)$$



which is brilliant as its a single closed-form, its differentiable across the whole domain and it is fairly easy to differentiate. The corr. derivative for the logistics function can be expressed in a self-related form as $\frac{d}{dx}S(x) = S(x)(1-S(x))$ which is handy to express the gradient in a compact form

$$\Pi'(x) = \frac{d}{dx} \left [  S \left ( \frac{x}{\sigma} \right ) - S \left ( \frac{x-l}{\sigma}\right ) \right] =  S \left ( \frac{x}{\sigma} \right ) \cdot \left [ 1-S \left ( \frac{x}{\sigma} \right ) \right ] -  S \left ( \frac{x-l}{\sigma} \right ) \cdot \left [ 1-S \left ( \frac{x-l}{\sigma} \right ) \right ] .$$

<br>

<p align=center><img src="inductance_profile.png"></p>

*The upper plot shows the perturbation profile and the corr. gradient which peaks exactly at the coil edges. Main perturbations occur when the particle enters and exits the coil. The characteristic width $\sigma$ should be chosen $\sim1/5$ of the projectile length - $\sigma$ can be seen as a tuning factor, the higher this length is the smaller the perturbation will be and vice versa. Small characteristic lengths lead to higher slopes and more isntantaneous changes causing higher induction voltage peaks in the inductor.*

The final model for the solenoid inductance is then given by

$$ \mathcal{L}(x) = \mathcal{L}_\text{min} + \lambda \cdot (\mathcal{L}_\text{max} - \mathcal{L}_\text{min}) \cdot \Pi(x) $$

$$ \mathcal{L}'(x) = \lambda \cdot (\mathcal{L}_\text{max} - \mathcal{L}_\text{min}) \cdot \Pi'(x) $$

In the following it will be assumed that the inductance and its derivative are known at any time $\mathcal{L}_t, \mathcal{L}'_t$ as they can be computed in every time step from the instantaneous variables $I_t, x_t, v_t$ etc.

## 3.2. Theoretical Inductance Response

An alternative derivation for the inductance profile can be obtained from the elegant force comparison i. e. the Biot-Savart law and reluctance force definition from section 2.2. We can compare both absolute (1-dim.) force expressions to find an elegant way to derive the theoretical inductance slope 

$$

\Rightarrow F_\text{Reluctance} \overset{!}{=} F_\text{Biot-Savart} \newline  \newline

\Leftrightarrow \frac{1}{2} \mathcal{L}'(x) \cdot I^2 = F_\text{BS}

\newline s
\newline
$$

Finally we can solve for $\mathcal{L}'(x)$ and re-use the force which was pre-calculated with the Biot-Savart approach

$$
\boxed {\Leftrightarrow \mathcal{L}'(x) = \frac{2}{I^2} \cdot F_\text{Biot-Savart} (x; \gamma, \mu, m, \rho, N, l, I, R) }
$$


For the inductance profile simply integrate over the derivative which leads back to the cumulative work applied to the projectile


$$\Rightarrow \mathcal{L}(x) = \int_{-\infty}^{x} d\tilde{x} \mathcal{L}'(\tilde{x}) = \frac{2}{I^2} \int_{-\infty}^{x} d\tilde{x} F_{BS}(\tilde{x}) = \frac{2W_\text{mag}(x)}{I^2} + \mathcal{L}_\text{air}$$

the obligatory constant $\mathcal{L}_\text{air}$ which arises from the integration can be interpreted as the physical boundary when the projectile is far away i.e. $|x|\rightarrow\infty$ while $\frac{2W_\text{mag}}{I^2}$ only describes the total work to magnetize and accelerate the projectile. 

$$\Rightarrow \mathcal{L}(x) = \frac{2W_\text{mag}(x)}{I^2} + \mathcal{L}_\text{air} = \frac{2}{I^2} \cdot \vec{\mathcal{M}}(x) \cdot \vec{B}(x) + \mathcal{L}_\text{air} = \frac{2 \mu V}{I^2} \cdot \vec{H}^2_\text{sol}(x) + \mathcal{L}_\text{air}$$

Finally insert the derived solution for $\vec{H}$ from section 1.2 to find the expression for the coil inductance in dependence of the projectile position

$$\boxed{\Leftrightarrow \mathcal{L}(x) = \frac{\mu m N^2}{2 \rho l^2} \left [ \frac{l - x}{\sqrt{R^2+(l-x)^2}} + \frac{x}{\sqrt{R^2+x^2}} \right ]^2 + \mathcal{L}_\text{air}} $$

Note that this inductance expression requires no field length tuning, unlike the previous method.

## 3.3. Back-EMF Effects

Kirchhoff's second law ensures that energy is conserved as the source voltage is working against the total resistance and induced voltages i. e. $V_{cc} = I\mathcal{R} + V_L$. Despite the coil not having electrical resistance (or negligible) it can still become a significant load since the circuit must supply energy to push current against the induced back EMF. Thus the minus sign from Lenz's law can be absorbed into $V_L$ and we just have to care about the magnitude.

$$ V_L = \mathcal{L} \dot{I} + \mathcal{L}'Iv $$
$$ \Leftrightarrow V_{cc} - \mathcal{R}I = \mathcal{L} \dot{I} + \mathcal{L}'Iv $$
$$ \Rightarrow \dot{I} = \frac{V_{cc} - (\mathcal{R} +\mathcal{L}'v)I}{\mathcal{L}} $$

In the last step we separated the variables to show how the time derivative of $I$ is influenced by itself. In the stationary case, where $x$ is fixed, the upper identity breaks down to a non-homogeneous first-order differential equation $\dot{I} = C_1 I + C_0$, where $C_i$ are constant. This would only apply if we would measure the coil current while the projectile stands still or is far away - which is unrealistic in the presence of the coil's field. Since  $\mathcal{L}'$ and $\mathcal{L}$ depend on $x$, which transversely varies with $t$, we cannot use simple LDE methods for the non-linear dependence. Instead, we rely on numerical integration methods to compute the current.

For discrete timesteps $\Delta t$ the time derivative can be separated into the deltas

$$ \Rightarrow \Delta I = \frac{V_{cc} - (\mathcal{R} +\mathcal{L}'v)I}{\mathcal{L}}  \Delta t$$


This update happens in each quasi-stationary integrator step, where $I, \mathcal{L}'$ and $\mathcal{L}$ are constant within each $\Delta t$ interval and $I$ refers to the most recent observation i. e. $I_\text{old}.$ The homogeneous part of the equation indicates an exponential solution for $I$, and it shows how the change in current is dampened by the inductance in the denominator. Also, we find that the motional EMF leads to damping perturbation which scale linearl with the current, the spatial inductance response slope, and the projectile velocity.

With Kirchoff's loop rule (for components in series) and Ohm's law we can express the total resistance $\mathcal{R}$ beforehand

$$ \mathcal{R} = \mathcal{R}_\text{coil} + \mathcal{R}_\text{external} = \xi \cdot \frac{l_w}{A_w} + \mathcal{R}_\text{ext} = \xi \cdot \frac{8 N R}{\theta^2} + \mathcal{R}_\text{ext}$$

where $\xi$ is the resistivity of the coil wire, $l_w$ is the total length of the wire, $A_w$ the wire cross sectional area, $\theta$ is the wire thickness and $\mathcal{R}_\text{ext}$ is the optional external resistance of the circuit which bounds the current from above i. e. $I \lt V_{cc}  /\mathcal{R}$. $\mathcal{R}$ can be used throughout the whole integration process as it's constant.

### Capacitive Voltage Source

The best way to maximize the coil's magnetic field and it's reluctant force is via large currents as we showed $F \propto I^2$. A common practice is to use large capacitors as they are able to release large amounts of stored energy in very short time and can withstand high voltages. The source voltage varies as

$$ V_{cc} (t) = \frac{Q_{cc}(t)}{C} = \frac{Q_{0}-\int_{t_0}^{t}I(t)\cdot dt}{C} = V_{0}- \frac{\int_{t_0}^{t}I(t)\cdot dt}{C}$$

Where $Q_{0}$ is the initial charge of the cap and $V_{0}$ the corresponding charge voltage.
Again we have to integrate in a discrete manner 

$$\Rightarrow V_{cc} (t) \approx V_{0} - \frac{1}{C} \sum_{t}^{} I_t\cdot \Delta t $$

Numerically, we initialize the voltage $V_{t} \leftarrow V_0$ and after every integration iteration we subtract the small voltage changes based on the outflow/inflow of stationary charge increments

$$\boxed {V_{cc} (t) \approx V_t = V_{t-1} - \frac{I_{t-1}\cdot \Delta t}{C}}  \wedge  V_{t=0} \equiv V_0  $$

## 4. Kinematics

To accelerate the projectile along different stages we need to prevent magnetic breaking which happens since the direction of force is proportional to $x$. If the projectile passes the coils center at $x=0$ the sign flips and so does the force, which would lead to a deceleration. Consider the stage time $t$ of a single coil stage. The experiment starts at $t=0$ when the first stage is switched on. Each coil needs to be switched off at a specific latency $\tau_i \gt 0 : x_i(t_i=\tau_i) = 0$ where the index $i$ labels the current stage. To compute the optimal latencies we need to know the time-resolved trajectory by integrating the equations of motion.

An optimal integrator choice which minimizes the bias (to first order) and conserves energies even after large number of time steps is the *leap frog* integrator. We initialize the acceleration, position and velocity.

$$x_t=x_0;  a=a_{mag}(x_t)=\frac{F}{m}(x_t;  \chi, \mu, \rho, N, L, I, R);  v_{1/2}=v_0 + a\cdot \frac{\Delta t}{2} $$

where we define the start point (negative distance to coil), the start velocity (at first stage 0) and the acceleration as the force term evaluated at 0. This initialization step can be applied to any stage with the respective $x_0$ and $v_0$. The $v_0$ can be derived from the final velocity (at time $T$ - stage time) of the previous stage i. e. $v_{0,S} = v_{T,S-1}$.

The integration logic is to use the velocity vector from $t+\Delta t/2$ and use it as the corrected velocity vector from the current point at time $t$. This results in the optimal mean velocity approximation over the time interval $\Delta t$, provided that the interval is small enough - leading to an unbiased estimate in first order. The loop logic is

1. Update poistion (drift)
   
$$x_{t+1} \leftarrow x_{t} + v_{t+1/2} \cdot \Delta t$$

2. Update velocity (boost)

$$v_{t+3/2} \leftarrow v_{t+1/2} + a(x_{t+1})\cdot \Delta t$$

3. Stop condition

$$ \text{stop when } x_t = \frac{l}{2} \vee  x_t < \frac{l}{2} \wedge x_{t+1} > \frac{l}{2} \Rightarrow \tau = t$$





## Integration Algorithm For Projectile Trajectory

Below is an optimized integrator algorithm to compute the time-resolved dynamics for custom solenoid- and projectile properties, voltage source and variable initial conditions.

### 0) Initialization
---
Initialize the following values

- initialize the time parameter $t \leftarrow 0$
  
- start position $\Leftrightarrow x_t \leftarrow x_0$

- start velocity $\Leftrightarrow v_t \leftarrow v_0$

- start velocity $\Leftrightarrow a_t \leftarrow 0$

- current $\Leftrightarrow I_t \leftarrow 0$

Compute the "empty" coil inductance i. e. $\mathcal{L}_\text{air} =  \frac{ \mu \pi N^2 R^2 }{ l }$

Determine the resistivity $\mathbf{\xi}$ of the coil wire material (for copper $\sim 1.72\times10^{-8}\Omega m$).

Determine the demagnetizing factor $\gamma$ depending on the projectile shape.

Determine the wire thickness $\theta$ in mm, drawn from wire gauge.

Compute the total resistance $\mathcal{R} = \mathcal{R}_\text{sol} + \mathcal{R}_\text{ext} = 8\xi N R / \theta^2 + \mathcal{R}_\text{ext}$





### 1) Determining Coil Inductance and Slope
---
Calculate the momentary coil inductance 

$$ \mathcal{L}_t = \frac{\mu m N^2}{2 \rho l^2} \left [ \frac{l - x_t}{\sqrt{R^2+(l-x_t)^2}} + \frac{x_t}{\sqrt{R^2+x_t^2}} \right ]^2 + \mathcal{L}_\text{air} $$

and gradient

$$ \mathcal{L}'_t = 0 \text{ if } t=0$$
$$ \mathcal{L}'_t = \frac{2m}{I^2} \cdot a_{t-1} \text{ if } t>0$$



### 2) Compute Coil Current
---
From the inductance and it's slope compute the change in current and update the coil current

$$ \Delta I = \frac{V_{cc} - (\mathcal{R} +\mathcal{L}_t'v) \cdot I_{t-1}}{\mathcal{L}_t}  \Delta t , \text{ where }  I_{t-1} = I_0 \text{ for the first run}$$

$$ I_t \leftarrow I_{t-1} + \Delta I $$

If the voltage source is a pre-charged capacitor with capacitance $C$, update the voltage source accordingly

$$ V_{cc} \leftarrow V_{cc} - \frac{I_{t-1}\cdot \Delta t}{C}  $$

otherwise just proceed to the next step.

### 3) Compute leapfrog step for projectile
---
Knowing the exact current in the coil allows to finally compute the acceleration

$$a_t  \leftarrow  \frac{1}{m} \cdot F(x_t;  \gamma, \mu, \rho, N,L,I_t,R)$$

the result is used to compute the change in velocity (boosting)

$$\text{if $t=0$: }    v_{t+1/2} \leftarrow v_{0} + a_t\cdot \Delta t / 2$$
$$\text{else: }    v_{t+1/2} \leftarrow v_{t-1/2} + a_t\cdot \Delta t$$

and analogously update the position for the next iteration (drifting)

$$x_{t+1} \leftarrow  x_{t} + v_{t+1/2} \cdot \Delta t $$


### 4) Stop Condition
---
There are two types of stop conditions:

1. Stop when the projectile reaches the center of the coil:
$$\text{stop when } x_t = \frac{L}{2} \vee  x_t < \frac{L}{2} \wedge x_{t+1} > \frac{L}{2}$$

2. Stop at a hard coded time $T$:
$$\text{stop when } t \ge T$$

### 5) Loop
---
Increment the time $t \leftarrow t + \Delta t$ and go to 1)