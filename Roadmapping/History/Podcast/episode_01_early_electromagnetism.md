---
episode: 01
title: "Early Electromagnetism 1800-1860"
era: "1800-1860"
chapter: 01_early_electromagnetism_1800_1860
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 22
animations_cued:
  - hist_faraday_induction
status: draft
---

# Episode 1 — Early Electromagnetism (1800–1860)

> Companion dialogue script for [[01_early_electromagnetism_1800_1860]]. Same primary research, conversational form.

## Cold open

**Historian:** On 20 March 1800, a letter arrives at the Royal Society in London. It's from a chemist in Como, in northern Italy — Alessandro Volta — and it's addressed to Sir Joseph Banks. The letter describes, in great detail, a column of brass and zinc disks separated by cardboard soaked in salt water. Touch the top of the column with one hand and the bottom with the other, and you receive a mild but *continuous* electric shock. Not the spark of a Leyden jar. Not the brief discharge of a frictional machine. A current that does not stop.

**Experimentalist:** And to anybody who'd been working with static electricity for the previous fifty years, that's the strange part. Everything they knew was *transient*. A Leyden jar charges up, you discharge it, that's it. Volta's pile gives them something they'd never had before — a *source*.

**Historian:** Volta calls it the *pile* because that's what it looks like — a stack. Within months replicas are running across Europe. Within a year Nicholson and Carlisle in London have used one to decompose water into hydrogen and oxygen. Within seven years Humphry Davy has isolated sodium and potassium with one.

**Physicist:** And before we go any further, I should say what this whole series is about — including this episode. We're walking through the history of physics from 1800 to about 1965, and we're doing it twice. Once in the standard way you'd find in any textbook, and once through a particular mathematical framework — Tepper Gill's *proper-time* or *dual* framework. The thing to understand from the outset: we are not proposing new physics. The dual framework reproduces every experimentally confirmed prediction of standard relativity and standard quantum mechanics, to within current measurement precision. What we're doing is exploring how the same experimental record can be rewritten in a different mathematical convention, and whether that rewriting makes any of the derivations cleaner or any of the conceptual difficulties dissolve. Where the dual framing makes a different numerical prediction, we'll always say what precision regime could distinguish it from the standard prediction. Most of the time, that regime is far beyond what current experiments can reach.

**Experimentalist:** Which for *this* episode means: the comparison is null. Nothing in the period 1800 to 1860 distinguishes the two framings. Currents are slow, sources are at rest in the lab, and the corrections you'd compute under the dual theory are roughly one part in `10^15` or smaller. The galvanometers of 1832 cannot see that.

**Physicist:** Right. So this episode is mostly establishing the historical record — what was discovered, who discovered it, what they had wrong, what they had right. The genuine divergence between standard and proper-time physics doesn't show up until Maxwell's 1865 paper and Michelson–Morley in 1887. Those are next episode.

**Historian:** Good. Then let's begin in Italy in 1800.

## Historical sweep

**Historian:** Volta's pile is the *operational birth* of current electricity. But for twenty years after 1800, nobody connects it to magnetism. The magnetic compass had been used in navigation for eight centuries. Lodestone had been studied by William Gilbert in 1600. Static electricity and magnetism were known to be separate things — they pulled on different objects, in different ways, with different forces. Nobody held a current-carrying wire next to a compass and looked.

**Experimentalist:** *Nobody*. For twenty years.

**Historian:** Until April of 1820 — possibly during a lecture demonstration, the dating is disputed — Hans Christian Ørsted, in Copenhagen, holds a current-carrying wire over a compass needle. The needle swings — *perpendicular* to the wire. Not parallel, which is what most people would have expected from a "magnetic" effect. Perpendicular. Ørsted writes it up in Latin, four pages, and distributes copies of [[orsted1820_experimenta_circa]] across Europe in July.

**Experimentalist:** A four-page paper. Latin. Self-published. And it triggers an avalanche.

**Historian:** Within weeks. The pamphlet reaches Paris in September. André-Marie Ampère reads it on the eleventh. By the eighteenth, exactly one week later, he is at the lectern of the Académie des Sciences presenting his first results. Two parallel current-carrying wires attract when their currents flow in the same direction. Antiparallel currents repel. Reading continues on the twenty-fifth of September and the second of October. The work becomes [[ampere1820_action_mutuelle]].

**Physicist:** And the *form* of the force law, the way Ampère writes it — force per unit length scales as the product of the two currents divided by the distance between the wires. He works it out for finite-length wire elements as a double integral, with the modern form of a $\,d\boldsymbol{\ell}_1 \times (d\boldsymbol{\ell}_2 \times \hat{r})$ vector triple product. It's the *Lorentz force on a current*, written down forty years before the Lorentz force itself.

**Historian:** And in the same weeks, Jean-Baptiste Biot and Félix Savart in Paris are measuring the *field* — not the force, the field — of a single long straight current-carrying wire. They find it scales as one over the distance from the wire. Two papers, [[biot_savart1820_magnetisme_pile]] and a 1821 follow-up, give us what we now call the Biot–Savart law. Ampère describes the force a current feels; Biot–Savart describes the field a current produces.

**Experimentalist:** Two halves of the same picture. And both done in three months.

**Historian:** Ampère also proposes — this is conceptually decisive — that *magnetism itself* is just the manifestation of microscopic, persistent molecular currents. The bar magnet is full of tiny invisible current loops. He calls them *Ampèrian currents*. It's the first attempt to reduce one phenomenon to the other. It survives, in a heavily modified quantum-mechanical form, all the way to today.

**Physicist:** Quick reality check from the proper-time side: everything we just described — Ørsted, Ampère, Biot–Savart — is steady-state magnetostatics from currents at rest in the lab frame. The collaborative speed $b = \sqrt{c^2 + u^2}$ in this regime is, to any practical precision, just $c$. The two framings give the same force laws and the same field magnitudes. This is the null-comparison era.

**Experimentalist:** Tagged "Gill-silent" in the chapter notes.

**Historian:** All right — jump forward seven years. Georg Simon Ohm, working in Bavaria with what he describes as woefully inadequate equipment, publishes [[ohm1827_galvanische_kette]]. It's a book — *Die galvanische Kette, mathematisch bearbeitet*, "the galvanic circuit, mathematically treated". And it derives, and experimentally verifies, what we now call Ohm's law: the current through a metallic conductor is proportional to the voltage across it.

**Physicist:** $V = IR$.

**Historian:** $V = IR$. Ohm constructs the analogy deliberately — he's been reading Fourier's *Théorie analytique de la chaleur*. Heat flux is proportional to the temperature gradient. So electric current, he reasons, must be proportional to the *potential* gradient. He's right. And the constant of proportionality is a property of the conductor alone — the *resistance*.

**Experimentalist:** And the reception?

**Historian:** Disastrous. The German physics establishment finds the work "speculative". Ohm is forced to resign his teaching post. He spends the next decade in academic exile. Recognition arrives slowly — first in France and Britain, eventually from the Royal Society, which awards him the Copley Medal in 1841. By then the telegraph industry is buying his book by the thousands of copies because nobody else has written down the rules for designing a long-distance signal line.

**Physicist:** Eponymy is the ultimate vindication. The unit of resistance is named after him.

**Historian:** Now we come to the central experimental discovery of the era. Michael Faraday, at the Royal Institution in London, has been brooding for ten years over a question. Ørsted showed that electricity *produces* magnetism. The converse must surely be true: magnetism must produce electricity. Faraday tries, and tries, and tries — for ten years he tries. Steady magnetic fields produce nothing. Permanent magnets next to wires produce nothing.

**Experimentalist:** And then on August 29th, 1831, in his laboratory beneath the Royal Institution, he wraps two coils of wire on opposite sides of a soft-iron ring. He attaches a galvanometer to one coil and a battery to the other. He closes the battery circuit. The galvanometer needle *kicks*. Once. He opens the battery circuit. The needle kicks back. With the battery connected and current flowing steadily, the needle reads zero.

**Historian:** Only *change* induces a current. That's the discovery. The published paper, [[faraday1832_experimental_researches_i]], appears in *Philosophical Transactions* in March 1832, and it's the first of an extraordinary thirty-volume series called *Experimental Researches in Electricity* that will run for the next thirty years. Section ninety of that first paper introduces the "Faraday disc" — a copper disc spinning between the poles of a permanent magnet, with brushes at the centre and the rim. The disc generates a continuous DC voltage proportional to its rotation rate. The first electromagnetic *generator*.

**Physicist:** And here's where I want to pause for the physics. We now write Faraday's discovery as $\mathcal{E} = -d\Phi/dt$. The EMF around a closed loop equals the negative of the rate of change of magnetic flux through it. Faraday himself never writes it in that form — the mathematical formalisation has to wait for Maxwell in the late 1850s — but everything we need to derive it is in Faraday's prose. The integral is implicit in his loop-and-flux language. The minus sign comes from a separate paper we'll get to in a moment.

**Experimentalist:** And independently — independently! — Joseph Henry at the Albany Academy in upstate New York is making the same discovery during the same summer. Probably slightly *earlier*, by a few weeks. But Henry has teaching obligations and doesn't publish until the following year, which is [[henry1832_currents_and_sparks]]. The priority goes to Faraday on publication date, and Henry ungrudgingly acknowledges it.

**Physicist:** Henry's distinctive contribution is *self-inductance* — that "kick" you feel, that spark you see, when you break the circuit of a long-wound coil. The coil's own changing field induces an EMF back in itself, opposing the change. The SI unit of inductance is the henry. So Faraday's work and Henry's work are complementary: Faraday on mutual induction, Henry on self-induction.

**Historian:** And the sign? Two years after Faraday, in 1834, Heinrich Lenz in St Petersburg publishes [[lenz1834_richtung_strome]] — "On the determination of the direction of galvanic currents excited by electrodynamic distribution". Lenz's law: the induced current is always in the direction such that its *own* magnetic field opposes the change that induced it.

**Experimentalist:** Which is conceptually decisive. With the *opposite* sign, you'd have a feedback loop that amplifies its own input. Energy conservation fails. You'd build a perpetual-motion machine.

**Physicist:** Right. Lenz's law is, in 1834, effectively a statement of energy conservation applied to electromagnetism — thirteen years before Helmholtz writes down the general principle in *Über die Erhaltung der Kraft*.

**Historian:** Around the same time, at the University of Göttingen, Carl Friedrich Gauss is doing something quite different but ultimately just as important. His paper [[gauss1832_intensitas]] introduces the idea of *absolute units* for magnetic field strength — expressing the Earth's magnetic field in centimeters, grams, and seconds alone, rather than in some instrument-specific calibration. He and Wilhelm Weber found the *Magnetischer Verein* — the Magnetic Union — in 1834, which coordinates simultaneous magnetic-field observations across observatories all over Europe using those absolute units.

**Physicist:** Which is the conceptual ancestor of the modern SI system. The methodological move — *expressing electromagnetic quantities in pure mechanical units* — is what makes nineteenth-century electromagnetism *quantitative* rather than instrument-specific.

**Experimentalist:** It's also what eventually leads to the ratio-of-units measurement that gives Maxwell the speed of light. Wilhelm Weber and Rudolf Kohlrausch will measure that ratio in 1856 and find it numerically very close to $3 \times 10^{10}$ cm/s.

**Historian:** And in 1857, Gustav Kirchhoff in Heidelberg uses that ratio to derive something extraordinary. [[kirchhoff1857_bewegung_elektrizitat]] — Kirchhoff's analysis of telegraph-cable signal propagation — shows that an electric disturbance travels along a wire at very nearly the speed of light in vacuum. Not because Kirchhoff is thinking about optics. He is solving a telegraphy problem. The match is left as a striking empirical coincidence.

**Physicist:** Kirchhoff does not draw the optical conclusion. That is reserved for Maxwell, eight years later. But the empirical anchor is here, in 1857. The number $c$ appears in an electromagnetic equation for the first time.

## Physics deep dive

**Physicist:** Let me take us into one derivation in detail — Faraday's induction law, the way Maxwell will eventually write it. Set up a closed conducting loop $C$ enclosing a surface $S$. The magnetic flux through that surface is $\Phi = \int_S \mathbf{B} \cdot d\mathbf{A}$. Faraday's experimental observation is that the EMF induced around the loop — the line integral of the electric field around $C$ — equals the negative time derivative of the flux:

$$\oint_C \mathbf{E} \cdot d\boldsymbol{\ell} = -\frac{d\Phi}{dt}$$

`[cue: animation hist_faraday_induction]`

**Physicist:** The animation we're cuing here shows the visual content: field lines threading through a translating loop; the loop moves into a region of stronger field; the flux through the loop increases; the induced EMF drives a current whose own magnetic moment *opposes* the increase. That last bit is the minus sign — that's Lenz.

**Experimentalist:** And practically — what does Faraday see in his actual ring experiment? He has two coils on a soft iron ring. Closing the primary circuit causes the field threading the secondary to jump from zero to some value $B_0$ very quickly. The secondary's flux changes from zero to $\Phi = B_0 A$. The induced EMF integrated over that brief interval is some finite total impulse. The galvanometer is ballistic — it doesn't show the EMF, it shows the *integrated charge* deflection. So you see a single quick kick.

**Physicist:** Right. And steady-state operation gives no kick at all, because $d\Phi/dt = 0$. That's the critical bit — the bit Faraday spent ten years missing. He kept trying with static, steady configurations, because that's how you'd think about electrostatics. He had to set up a *changing* configuration before he saw anything.

**Experimentalist:** And the dynamo — the rotating disc — does that, mechanically and continuously. Mechanical motion through a non-uniform field gives a continuous $d\Phi/dt$, which gives a continuous EMF, which drives a continuous DC current.

**Physicist:** Which is the prototype of every electric generator built since. The hydroelectric dam. The nuclear power station. The bicycle dynamo. All of them are extended Faraday-disc descendants.

## Proper-time interlude

**Physicist:** Now let me contrast the standard derivation with how the dual framework would frame this. In Tepper Gill's *Two Mathematically Equivalent Versions of Maxwell's Equations* — that's the paper [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] in our verification corpus — there's a substitution that takes the ordinary time derivative $\partial_t$ to a proper-time derivative $\partial_\tau$, weighted by the collaborative speed $b = \sqrt{c^2 + u^2}$. Specifically, $\partial_t = (b/c) \,\partial_\tau$.

**Experimentalist:** Wait — *whose* proper time? If the source is at rest in the lab frame, there's only one clock.

**Physicist:** That's exactly the point for this era. When everything's at rest in the lab — Faraday's coils don't move, Henry's coils don't move, only the field configurations change — then $u = 0$ and $b = c$ identically. The dual induction law is exactly $\mathcal{E} = -d\Phi/d\tau$, which under $u = 0$ is exactly $\mathcal{E} = -d\Phi/dt$. The two framings collapse to the same equation. Faraday's experimental record is reproduced trivially.

**Experimentalist:** So when *would* they diverge?

**Physicist:** They diverge when the source moves at a non-negligible fraction of $c$. The correction is roughly $1 + u^2/(2c^2)$. For the rotor tip of Faraday's 1831 disc — copper, hand-cranked, maybe ten metres per second — that's a correction of one part in `10^15`. The galvanometers of 1832 had a precision of maybe one part in `10^4`. The dual-theory correction is eleven orders of magnitude below the experimental floor.

**Experimentalist:** Which is what makes this era a null comparison.

**Physicist:** Exactly. Same for Kirchhoff's telegrapher equation. The signal velocity on a Victorian copper transmission line was a tiny fraction of $c$ — partly because the dielectric was glass or gutta-percha, partly because resistive losses dominated. The dual-theory correction to the dispersion relation is buried below the noise floor by something like ten orders of magnitude.

**Experimentalist:** So when do we see divergence?

**Physicist:** Next episode. Maxwell's 1865 *Dynamical Theory of the Electromagnetic Field* introduces the displacement current, which makes light propagate as an EM wave. Once light is in the picture, the relevant velocity *is* $c$. Now $b$ and $c$ are at the same scale. The Michelson–Morley apparatus of 1887 is designed to detect that scale. And — this is the part that makes the campaign worth doing — Gill's dual formulation reframes the failure of the aether to be detected without invoking length contraction. We'll get to that.

## Closing

**Historian:** So that's the period 1800 to 1860 in one episode. Volta to Kirchhoff. The pile, Ørsted, Ampère, Biot–Savart, Ohm, Faraday, Henry, Lenz, Gauss, Kirchhoff. Ten primary sources, a handful of retrospective ones, all of it in the bibliography for this chapter. The full reference list is in the show notes, with wikilinks into the chapter document for anyone who wants the citations directly.

**Experimentalist:** And the methodological note for listeners who came in expecting the dual-theory contrast: there isn't one this episode. Not yet. Sources don't move fast enough, instruments don't measure precisely enough. Both framings predict the same outcomes for every experiment in this era. The genuine divergence is in 1887, with Michelson and Morley, and we'll get to that.

**Physicist:** And the reminder we set up with at the start: this entire series is exploring mathematical conventions for the same experimental record. We're not proposing new physics. Where the dual framing predicts something different from standard relativity or standard quantum mechanics, we will always say what precision regime distinguishes them. Most of the time, that regime is far beyond current measurements. The cases where it *is* within reach are precisely the load-bearing payoffs of the verification campaign that motivates this whole project — and they sit in chapters two through five.

**Historian:** Next episode: Maxwell's 1855 derivation of Faraday's lines of force in mathematical form, the 1865 *Dynamical Theory* with the displacement current and the prediction $c = 1/\sqrt{\mu_0 \varepsilon_0}$, Hertz's 1888 confirmation that EM waves *exist*, and Michelson–Morley's null result the year before. The first chapter that's not a null comparison. Thanks for listening.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`; primary + retrospective sources cited above.
