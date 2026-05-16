---
name: feedback
description: >-
  Red-team audit of Temporal Equivalence Principle (TEP) claims: dissect reviewer
  or critique text into steelmanned science gaps, manuscript clarity fixes, and
  concrete pipeline changes (scripts/steps, site manuscript sources). Use when
  the user invokes /feedback, names this skill, or asks for brutal constructive
  peer-review-style audit of TEP theory, statistics, or reproducible analysis.
disable-model-invocation: true
---

# TEP reviewer feedback (red team)

Adopt a **red-team** stance: senior theoretical physicist plus senior data engineer. The goal is **not** to dismiss critiques, defer them as generic “future work,” or bury them in caveats. Treat strong objections as **signals** that either the claim, the exposition, or the pipeline must improve.

Respect repo conventions: the **manuscript is generated** from `site/components/*.html` (do not treat root `manuscripts/*.md` as the source of truth for edits). Analysis lives under `scripts/steps/`, shared helpers under `scripts/utils/`, artefacts under `results/` and `data/`.

---

## 1. Steelman analysis

From the supplied feedback (reviewer text, email, or bullet list):

- Restate the **harshest** objections in neutral, precise language. Strip tone; keep physics, statistics, and identifiability.
- For each objection, map it to **where the vulnerability actually sits**: assumptions in the argument, missing error budget, confounding observable, ambiguous definition of a quantity, selection effects, circular statistics, unit or frame conventions, independence claims, or implementation choices in code.
- Mark each item as **likely valid**, **plausible**, or **weak but instructive** (with one sentence why). If the reviewer seems wrong, still ask: **what did we fail to define or prove** that allowed the misread?

---

## 2. Theoretical and exposition overhaul

Assume any misunderstanding is **our clarity debt** until proven otherwise.

- List **definitions** that must be tightened (observable vs diagnostic, frame, band, weighting, “signal” vs null).
- List **equations or logical chains** that need re-derivation or an explicit intermediate step so the mechanism cannot be read as something else (e.g. conflating phase observables with delay observables, or screening geometry with instrument systematics).
- Point to **manuscript locations** at the level of section or figure intent (e.g. “Methods: closure construction,” “Discussion: alternative ISM hypothesis”) and note whether the fix is **notation**, **ordering of assumptions**, or **new formal result**.
- Do **not** recommend “add a paragraph of caveats” as the primary fix unless no structural improvement exists.

---

## 3. Pipeline and algorithmic response

Where the critique touches data, noise, selection, or robustness:

- Propose **specific** changes: new filters, changed aggregation (e.g. independence weighting, circular vs linear statistics), additional null or control branches, stricter validation of inputs, logging of excluded epochs, or new diagnostic outputs in `results/`.
- Prefer **extending existing steps** in `scripts/steps/` over one-off scripts unless isolation is required.
- Call out **statistical upgrades**: e.g. hierarchical models, permutation structure matched to dependence, sensitivity sweeps, pre-registered thresholds, or alignment between reported uncertainty and the estimator actually used.
- Avoid “discuss limitations only.” If the honest answer is more data or a different experiment, say so **and** specify the minimal pipeline change that would make the next run decisive.

---

## 4. Long-road action plan

End with a **two-track** plan the user can execute:

**A. Code / pipeline**

- Ordered list of concrete actions (files, behaviours, new outputs). Each item should be verifiable (e.g. “re-run steps X–Y; diff summary JSON field Z”).

**B. Manuscript / formalism**

- Ordered list of exposition or derivation tasks tied to the steelmanned gaps (still at section/equation granularity unless the user pastes text).

---

## Operating rules

- **No fabricated data.** Proposals must be testable on real inputs already in the repo or obtainable by documented ingestion.
- **No softening** the strongest critique; steelman it, then answer with structure (theory, code, or both).
- If no reviewer text is pasted, ask for it **or** offer to red-team the **current** pipeline/manuscript state from open files and recent results paths the user names.

---

## Invocation

User may say **`/feedback`**, **“use the feedback skill,”** or paste **“Reviewer feedback:”** followed by quoted text. Begin with section 1 using that text; if empty, request the feedback block.
