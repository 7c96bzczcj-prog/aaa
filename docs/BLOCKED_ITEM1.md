# Item 1 (the two numbers) — BLOCKED, and exactly what unblocks it

**Status: cannot be done in this environment.** Not deferred again — unreachable.

## What was checked

```
find / -iname "*dean2024*" -o -iname "*falsify*"     -> nothing
grep -rli "oxphos|seahorse|ccr5|MC38" (csv/md/py/tsv) -> only this repo's own T1 table
ls /home/user/                                        -> aaa/  (this repo only)
```

This container holds **one** repository, `7c96bzczcj-prog/aaa`, the NK
lineage-assertion screen. There is no `dean2024/`, no `falsify/`, no OXPHOS
analysis, no MC38 data.

## The two numbers, restated so they are not lost

1. **Unique values of the `donor` column in the OXPHOS analysis.** If those turn
   out to be library IDs (FCA-style) rather than genuine donors, the
   "11/11 donors, p = 0.001" result is void — pseudoreplication at the donor
   level, the same failure as `n = wells`. That result is one of three independent
   supports for the main line.
2. **Detection rate of `Ccr5` in MC38.** Open
   `dean2024/falsify/E_pos_frequency_vs_level.csv` first — the number may already
   be there and not need a re-run.

**Follow-on if number 1 comes back as library IDs:** grep that workspace for
archived p-values equal to `2/2^k` — those are the signature of a permutation test
whose null was built from too few distinct units, and they will cluster in the
same analyses.

## What unblocks it — any one of these

1. **Attach the workspace to this session** (if it is a GitHub repo, I can add it).
2. **Paste `E_pos_frequency_vs_level.csv`**, or just its header plus the `Ccr5`
   row — that alone answers number 2.
3. **Paste the `donor` column's unique values**, or the metadata table that
   defines it — that answers number 1 and tells us whether the grep is needed.

Option 3 is a few lines of text and resolves the item that gates the proposal
structure.
