# Endogenous Inflation in a Heterogeneous-Agent Model of Debt-Financed Fiscal Stimulus

*A note on why the extension exists, what it adds, and what it does and does not deliver.*

---

## 1. The problem with feeding inflation

The baseline model is an open-economy heterogeneous-agent economy with capital, long-term nominal
government debt, foreign demand for that debt, and HSV progressive labour taxation, in which
households have Epstein–Zin preferences over a GHH consumption–labour composite. It is used to study
the 2020–24 US fiscal expansion: a large, debt-financed transfer programme, financed at the margin by
issuing nominal government bonds, whose real burden depends on the inflation that follows.

In the original formulation the inflation path is **exogenous**. The observed sequence is fed in
quarter by quarter, and everything nominal follows from it. This is a defensible modelling choice for
the positive exercise — it lets the model match the observed debt revaluation exactly — but it has a
sharp consequence for the counterfactuals.

Because the price level is imposed, fiscal policy cannot move it. And because wages are flexible,
employment is whatever the household's intratemporal first-order condition delivers at the prevailing
real wage and retention rate: there is no channel from aggregate demand to the labour market. Fiscal
policy therefore acts on quantities only through the real interest rate and the tax schedule.

The bite comes in the counterfactual section. Each counterfactual changes the size or composition of
the fiscal package while holding the inflation path at its observed values. So a counterfactual in
which transfers are roughly tripled still has 9.06% annualised inflation in 2022Q2, and the real value
of the debt issued to finance those transfers is eroded by exactly as much as in the data. The
mechanism the paper is built around — that a debt-financed expansion is partly paid for by inflating
away nominal liabilities — is switched off precisely in the experiments designed to interrogate it.
Put plainly: **the counterfactuals inherit the observed inflation while changing the fiscal policy that
would have produced it.**

The extension makes inflation endogenous so that this objection can be answered rather than deflected.

---

## 2. What has to be added, and why it is exactly three things

It is worth being precise about the accounting, because the temptation is to add too much.

In the original model, each transition period contributes one market-clearing condition — asset market
clearing, equating household savings to the capital stock plus the domestically held government debt —
and one unknown, the real interest rate. The wage follows from the firm's first-order condition given
the rental rate, the labour tax retention level adjusts residually to satisfy the government budget,
and hours come from the household's first-order condition. The system is determined.

**Sticky wages add one equation and one unknown.** The new unknown is the degree to which hours are
demand-determined rather than chosen; the new equation is the wage Phillips curve that governs it.
Balance is preserved.

**The monetary block adds one net equation, not two.** This is the step that is easy to miscount. A
Taylor rule is one equation, but it also introduces the *nominal* interest rate, which did not appear
anywhere in the real model. The Fisher equation is a second equation, but its role is to eliminate that
new variable, tying the nominal rate to the real rate and expected inflation. Two equations, one new
auxiliary variable: **net one equation**, matching the one genuinely new unknown, the inflation rate.

So the extension is: hours become an unknown governed by a wage Phillips curve; inflation becomes an
unknown governed by a monetary rule closed through Fisher. Nothing else changes. In particular the
household problem, the firm problem, the government budget and the foreign block are untouched.

---

## 3. The union block

**In words first.** Labour is supplied to a continuum of unions, each of which organises workers into a
differentiated labour variety and posts a nominal wage on their behalf. Households supply whatever
hours are demanded at the posted wage. Unions face a quadratic cost of changing the nominal wage, so
they adjust it only partially in response to shocks — the standard Rotemberg device, here applied to
wages rather than prices.

A union setting its wage weighs two things. Raising the wage raises its members' earnings, valued at
the after-tax marginal rate that the progressive tax schedule delivers, but it reduces hours demanded.
Not raising it leaves members working hours whose marginal disutility may exceed what they are paid at
the margin. The gap between those two objects — the marginal rate of substitution between leisure and
consumption on one side, the after-tax marginal wage on the other — is the **labour wedge**.

The Phillips curve says: when the wedge is positive, so that members are working more than they would
willingly supply at the going after-tax wage, unions push nominal wages up, and they do so more
aggressively the larger the wedge and the cheaper it is to adjust. Wage inflation today also depends on
expected wage inflation tomorrow, discounted, because the adjustment cost is forward-looking.

**In symbols.** Writing `piw` for nominal wage inflation and `pi_bar` for the indexation target,

```
(piw_t − pi_bar)(1 + piw_t) = kappa_w · wedge_t + beta · (piw_{t+1} − pi_bar)(1 + piw_{t+1})
```

where the wedge aggregates, across the wealth and productivity distribution, the difference between the
marginal disutility of the hours actually worked and the after-tax marginal value of the earnings they
generate, each household weighted by its marginal valuation of the consumption composite. Price
inflation follows from wage inflation and the change in the real wage, since the real wage is nominal
wage over price level by definition.

**What `kappa_w` means.** It is the slope of the wage Phillips curve, inversely proportional to the
adjustment cost. A small `kappa_w` is a flat curve: wages barely respond to the wedge, so hours absorb
shocks and the labour market is far from its flexible-wage allocation. A large `kappa_w` is a steep
curve: wages move quickly, the wedge is driven towards zero, and hours return to the household's
first-order condition. In the limit `kappa_w → ∞` the wedge is exactly zero and the model collapses back
to the flexible-wage original. This limit matters later: it means the sticky-wage extension is a
one-parameter family that *contains* the original model as an endpoint, so any result that holds
uniformly across `kappa_w` is a statement about the whole family.

One modelling choice deserves flagging. The steady-state wage markup is assumed to be offset by a
constant employment subsidy. This is conventional, and it has a convenient consequence here: the
elasticity of substitution across labour varieties cancels out of the Phillips curve entirely, and the
adjustment-cost parameter appears only through `kappa_w`. The extension therefore introduces exactly
one new structural parameter.

---

## 4. Why hours are rationed proportionally

When wages are sticky, hours are demand-determined, and one must say *how* the shortfall in hours is
distributed across a heterogeneous population. The choice made here is **proportional rationing**: every
employed household works the same fraction of its steady-state hours,

```
n_it = s_t · ntilde_z
```

where `ntilde_z` is household type *z*'s steady-state hours and `s_t` is a single economy-wide scalar.
The alternative — uniform rationing, in which everyone loses the same *number* of hours — is at least as
natural a priori. Three considerations favour the proportional rule, and the third is decisive.

**It preserves the cross-sectional shape of earnings.** Under proportional rationing every household's
labour earnings are scaled by the same factor `s_t`, so the distribution of log earnings shifts but does
not change shape. The income process is estimated to match that cross-sectional distribution; leaving
its shape intact means the process does not have to be re-estimated for the sticky-wage economy, and
the calibration carries over unchanged. Uniform rationing compresses the earnings distribution — it
takes proportionally more from low-hours households — and would require re-estimation.

**It is exact in this model rather than approximate.** With proportional rationing, the flexible-wage
value of `s_t` — the value that would obtain with no wage rigidity — is a closed-form function of the
real wage and the retention level alone, and it is **independent of productivity type** to within
numerical precision (verified to seven parts in 10^16). A single scalar is therefore not an
approximation to a type-specific object; it is the exact solution.

**Under HSV progressivity it delivers a type-independent wedge, and uniform rationing would not.** This
is the substantive reason. With progressive taxation the marginal retention rate depends on the level of
earnings, so the after-tax marginal wage differs across households. Setting the household's first-order
condition in the sticky-wage economy against the same condition in steady state, the productivity level
and the steady-state hours profile enter *both* the marginal-disutility side and the marginal-retention
side in exactly the same way, and cancel. What survives is a relation between `s_t`, the real wage and
the retention level that holds identically for every type. Under uniform rationing the hours term no
longer scales with `ntilde_z`, the cancellation fails, and the required adjustment becomes
type-specific — so a single labour wedge, and hence a single Phillips curve, would no longer be
well defined without further approximation.

---

## 5. Why the exogenous sequence is a monetary shock, not a cost-push shock

Once inflation is endogenous, the model will not reproduce the observed inflation path on its own.
Something exogenous must be added to make the baseline match the data. There are two natural candidates
— a residual in the Phillips curve (a cost-push shock) or a residual in the monetary rule (a policy
shock) — and the choice is not innocuous.

**The counting argument.** A free shock in the Phillips curve makes that equation non-binding. Whatever
pair of hours and wage inflation the rest of the model produces, the residual can absorb the difference,
so the Phillips curve stops determining anything. The equation is added and then immediately
neutralised, leaving hours undetermined. A residual in the monetary rule does not have this problem: the
Phillips curve still ties hours to the wedge, and the residual only affects how the central bank maps
inflation into the nominal rate.

**The structural argument, which is the more important one.** Suppose the inflation path were imposed
and the rule intercept held fixed. Then the Taylor rule delivers the nominal rate mechanically from
observed inflation, and the Fisher equation delivers the real rate from the nominal rate and observed
expected inflation. The real interest rate is thus **pinned by the monetary block before asset market
clearing is ever consulted**. But asset market clearing must still hold. With the real rate already
determined, the burden of clearing the asset market falls on the other endogenous quantity in that
market — employment, and through it output and savings.

That inverts the economics of the model. In the paper's baseline the real interest rate rises during
the expansion *because* the government is absorbing savings and crowding out capital; that is the
content of the real-rate path in Figure 4(e), and it is what makes the welfare incidence fall on
different households differently. If the real rate is instead pinned by the monetary rule, that
crowding-out channel is gone: the real rate follows the policy rule and employment adjusts residually.
The mechanism the paper documents would be an artefact of the exogeneity assumption rather than a
result.

So the exogenous sequence enters the monetary rule. The real rate continues to be determined where it
should be, in the asset market, and the policy residual absorbs the discrepancy between the rule and
the inflation the data displays.

---

## 6. What the extension delivers, and what it does not

**It delivers an endogenous employment response.** With sticky wages, fiscal policy moves hours through
a channel absent from the original model. A change in the capital tax, for instance, does not enter
labour supply directly and so cannot move employment in the flexible-wage economy at all; under the
union block it moves the real interest rate, hence the real wage, hence the marginal retention rate,
hence the labour wedge, and employment responds. In the counterfactual experiments this response is on
the order of one to two percent of hours — small, but no longer identically zero, and it is the
right sign.

**It delivers an endogenous inflation response, and that response is small.** With the estimated policy
rule and no monetary shocks, the model generates inflation that never rises more than a few hundredths
of a percentage point above target, against an observed peak of seven percentage points above target.
Making the wage Phillips curve six times steeper roughly triples that number and leaves it still under
one percent of what the data show. Because the flexible-wage economy is the `kappa_w → ∞` endpoint of the
same family, this is not a statement about a particular calibration of wage rigidity: no degree of wage
stickiness within the family generates the surge.

**It does not attribute the 2021–22 inflation to fiscal policy, and one should be candid about why.**
Prices are flexible and firms are competitive, so there is no direct channel from aggregate demand to
the price level. The *only* nominal rigidity in the model is in wages. Inflation must therefore travel a
long chain: the fiscal expansion crowds out capital, which moves the real interest rate, which moves the
real wage through the firm's first-order condition, which moves the labour wedge, which moves wage
inflation, which finally moves price inflation. Every link in that chain attenuates, and the
quantitative result is that the real rate barely moves at all — the measured elasticity of the real rate
to inflation is of order a few thousandths — so almost nothing reaches the price level. Compounding
this, a Taylor rule that satisfies the Taylor principle actively offsets whatever pressure does arise;
with no monetary shocks, anchoring inflation at target is close to what such a rule is defined to do.

It is also worth stating what is simply absent. The model has no imported energy or goods prices, no
markup or cost-push disturbances, and no supply-chain frictions. These are the leading candidate
explanations for the 2021–22 episode in the empirical literature. The finding is therefore best read
narrowly: *this* model's single inflation channel — a wage Phillips curve driven by a labour wedge —
does not produce the surge. That is a weaker and more defensible claim than "fiscal policy did not cause
the inflation."

---

## 7. How to describe the policy residual

The exogenous monetary sequence should be described as **the deviation from the estimated Taylor rule
required to support the observed inflation path**. It is roughly −4.2 percentage points annualised at the
2022Q2 peak, and negative values mean policy is more accommodative than the rule prescribes.

It should **not** be described as a measurement of what the Federal Reserve did, and the reason is worth
being explicit about in the text, because a reader will otherwise draw the wrong inference from its size.

The model has perfect foresight. Households and firms know the entire inflation path from the outset. In
2022Q2, with inflation running around 9% and the model's real rate around 1%, the Fisher equation
requires a nominal interest rate near 10% for the intertemporal problem to be consistent. The Federal
Reserve's policy rate at that date was around 2%. That eight-percentage-point gap is not evidence that
policy was extraordinarily loose; it is a consequence of assuming agents foresaw the inflation that in
fact surprised them. Under perfect foresight, an inflation path can only be an equilibrium if nominal
rates rise to meet it.

The residual reported is the gap relative to the *rule*, not relative to the observed policy rate, and
it is a smaller number for that reason. But it inherits the same interpretive caveat. It is an
accounting device that quantifies how far the model's monetary block must depart from its estimated rule
to be consistent with the observed price path — a measure of what the model cannot explain, not an
estimate of a policy shock. Read that way it is informative: it says that essentially the whole of the
2021–22 inflation, in this model, sits outside the fiscal-monetary interaction the model describes.
