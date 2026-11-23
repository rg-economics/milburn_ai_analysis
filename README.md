# Is the Milburn Review AI proof? Youth NEET model

RG Economics, November 2026. Code and data behind the blog [Is the Milburn Review AI proof?](https://www.rg-economics.com/milburn-review)

**All results are RG Economics modelling, not official figures.**

## What it does
The Milburn Review projects that its reforms would cut the youth NEET rate (the share of 16-24 year olds not in education, employment or training) from 13.4% to 6.8% by 2036/37. Its projections have no explicit AI channel. This model runs those projections through the five Government Office for Science (GO-Science) "AI Scenarios 2030" and tests whether Milburn's policies still work in each one (the "AI-proofing test").

## Run it
```
Rscript run_all.R
```
or one step at a time:
```
Rscript R/00_extract_a06.R   # takes the 16-24 series out of ONS table A06
Rscript R/01_data.R          # builds all the model inputs from data/raw
Rscript R/03_run.R           # fits the model, runs the scenarios, Monte Carlo, policy tests and sensitivities (about 15-20 min)
Rscript R/04_charts.R        # makes the charts in outputs/charts
```
`R/02_model.R` holds the model itself and is loaded by `03_run.R`.

R packages: readxl, readr, dplyr, tidyr, purrr, stringr, ggplot2, scales, systemfonts. Everything runs from published data with a fixed random seed (set in `assumptions/global.csv`).

## Change the assumptions
Everything you might want to change is in `assumptions/`:
- `scenarios.csv`: how each GO-Science scenario is turned into numbers (central, low and high values, with page references)
- `global.csv`: assumptions that apply to every scenario, such as today's AI effect on job openings, how much hiring depends on openings, the policy split, the pass mark and the number of runs
- `scenario_targets.csv`: how big the slump is in *Transformation Economy* and *Take-Off*
- `milburn_paths.csv`: Milburn's baseline and reform paths (his data tables, Table 3)

## How the model works
1. **Four states.** Every 16-24 year old is in work, in education or training, NEET and unemployed, or NEET and inactive. The model runs quarterly from Apr-Jun 2026. New 16 year olds join and 25 year olds leave.
2. **Flows between states.** The ONS doesn't publish these for young people, so the model starts from the all-age rates (ONS X02, latest eight quarters). Four numbers are fitted so the model reproduces today's numbers of young people in each state (ONS NEET table and A06): multipliers on job loss, job finding and unemployed-to-inactive moves, and the rate at which people leave education. There are four numbers but only three independent targets, so a small penalty picks between near-equal solutions. The results are in `outputs/tables/fitted_flow_parameters.csv`.
3. **Milburn's baseline.** How easy it is to find a job is set quarter by quarter so the run with no AI change follows Milburn's 13.4% baseline. This match is by design, not a test of the model.
4. **Job finding.** Job finding rises with the number of job openings and falls as more people look for work (a standard Cobb-Douglas matching function, with `alpha` setting how strongly hires respond to openings).
5. **AI channel.** 104 occupation groups (SOC 2020 minor groups) are weighted by how many 16-24 year olds work in them (APS 2025, ONS ad hoc 3537). Each gets AI exposure scores from the ILO (Working Paper 140, 2025), linked to UK occupations through the ONS SOC 2020 coding index. The ILO sorts occupations into four exposure grades. Which grades count as AI replacing work (which cuts hiring) and which as AI helping with work (which can add hiring) is RG Economics' judgement, not an ILO measure. AI cuts job openings in exposed jobs first. Some of the lost hiring reappears in less exposed jobs after a delay. Job losses among existing workers come later and are smaller. Each scenario also has a wider effect on demand (`D_demand` in `scenarios.csv`).
6. **Today's AI effect.** A 20% cut in job openings for a fully exposed occupation at today's level of AI use (range 10-30%), informed by Stanford's 19% employment gap for 22-25 year olds and King's College London's 23.4% fall in job adverts for highly exposed roles. Hires fall by roughly half as much as openings.
7. **The slump.** *Transformation Economy* and *Take-Off* describe recession-level unemployment by 2030, which AI in young people's own jobs doesn't reach on its own. So they are run three ways:
   - `jobs_only`: AI acts only through the jobs young people do.
   - `slump_fade` (central): an extra fall in youth job openings starts in 2029 and is sized in every run so the youth NEET rate in 2030/31 rises by as much as in 2008-11 (1.7 points for *Transformation Economy*, the average 2009-11 rise; 3.0 points for *Take-Off*, the rise to the 2011 peak). It then lifts in a straight line by early 2037.
   - `slump_held` (upper bound): as above, but the slump stays in place to 2036.

   Runs that would need openings to fall by more than 85% are dropped (131 of 2,000 *Take-Off* runs). Sizing the slump to the 2009-11 levels (15.6% and 16.9%) instead is a sensitivity (`sensitivity_targets.csv`).
8. **Milburn's policies.** Modelled as fixed numbers of places and hires each quarter, which is how his annex costs them: education and training places lasting about a year (`edp_spell_q`), subsidised hires whose take-up falls when openings fall (`beta_takeup`), and activation, health and welfare support whose job entries depend on job finding. A smooth profile is fitted so that, with no AI change, the package follows Milburn's reform path in every year (`reform_path_check.csv`). The education share (33%) is fitted so the package puts 440,000 young people into work by 2036/37. The rest is split 30:70 between subsidised hires and support, in line with their costs. Both matches are by design.
9. **Uncertainty.** 2,000 Monte Carlo runs per scenario, each drawing every assumption at random from its range. The bands show the spread of plausible assumptions, not probabilities. There is also a one-at-a-time sensitivity check, rank correlations and a breaking-point test.
10. **AI-proofing test.** A policy passes if it keeps at least 80% of the fall in the NEET rate it achieves with no AI change, in every scenario, at both the median run and the 10th percentile run.

## Outputs
- `outputs/tables/headline_2036_37.csv`: NEET rates in 2036/37 by scenario, with and without Milburn's policies
- `outputs/tables/policy_robustness.csv` and `ai_proofing_test.csv`: the AI-proofing test
- `outputs/tables/slump_size.csv` and `central_openings_and_slump.csv`: how big the slump is
- `outputs/tables/package_employment_effect_k.csv`: young people the package gets into work (thousands)
- `outputs/tables/policy_breaking_point.csv`: how big a fall in openings each policy can cope with
- `outputs/tables/sensitivity_*.csv`: sensitivity checks
- `outputs/tables/uk_aps_exposure_check.csv`: what the UK data shows so far
- `outputs/charts/`: the charts in the blog (chart 5 is shown in the blog for *Transformation Economy* with AI in jobs only; the other three versions are here too)

`mc_draws.rds` (every Monte Carlo run) is not included because of its size. It is recreated when you run the code.

## Main limitations
1. **Youth flow rates are fitted, not observed.** The fitted multipliers fill gaps in the data rather than describe real behaviour. LFS longitudinal data would fix this.
2. **Turning the scenarios into numbers is a judgement call.** GO-Science gives scores and descriptions, not model inputs, and stops in 2030. Everything after 2030 is RG Economics' assumption.
3. **The slump is the biggest single assumption** for *Transformation Economy* and *Take-Off*. It hits all youth jobs equally. It stands in for older workers who lose their jobs to AI competing for young people's jobs, which isn't modelled directly.
4. **Exposure scores are global**, not UK-specific. The SOC-to-ISCO link uses job title counts as weights. About 1% of youth employment uses filled-in (imputed) scores.
5. **The policy split is assumed.** Milburn doesn't publish NEET effects by policy, and his "system transformation effect" isn't modelled separately (see `sensitivity_lever_split_and_spell.csv`).
6. **Things the model leaves out.** More young people staying on in education during a downturn; subsidised hires replacing hires firms would have made anyway; long-term scarring; Milburn's Guarantees growing in a slump (they are entitlements but are modelled as fixed numbers); job quality and pay; population change (headcounts use the 2026 youth population of 7.54 million).
7. **Recovery is too slow.** A 2008-style shock raises the youth NEET rate by about as much as in 2008-11, but the model recovers much more slowly than the UK did after 2011.

## Data
See `DATA_LICENCES.md`. All raw data is public and included in `data/raw`.

## Licence
The code is released under the MIT licence (see `LICENSE`). You can use, change and share it as long as you keep the copyright notice. The data in `data/raw` keeps its own licences (see `DATA_LICENCES.md`).

Please cite as: RG Economics (2026), *Is the Milburn Review AI proof?*, https://www.rg-economics.com/milburn-review
