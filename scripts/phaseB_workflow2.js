export const meta = {
  name: 'nk-phaseB-decidability',
  description: 'Phase B: four decidability gates for every claim whose strongest Phase A evidence is INF_*',
  phases: [
    { title: 'Decide', detail: 'B1 lineage-decidable, B2 both-answers-actionable, B3 compartments' },
    { title: 'Occupancy', detail: 'B4 has anyone already done it with lineage/clonal methods' },
  ],
}

// args: [{claim_id, claim_text, evidence_type_strongest, load_bearing, notes}, ...]
const CLAIMS = args

const PREAMBLE = `You are running Phase B of a PREREGISTERED screen of NK-cell lineage claims.
Phase A already graded the evidence. Your job is ONLY decidability. Do not re-grade evidence and
do NOT offer any opinion on whether the question is scientifically valuable or worth doing -
that judgement is explicitly out of scope for this round.

The single primitive clonal barcoding provides: FOR ANY TWO CELLS, WHETHER THEY SHARE AN
ANCESTOR AND HOW RECENTLY. Nothing else. It does not measure function, phenotype causation,
or mechanism.

B1 lineage_decidable
  - "the progeny of A are B"                -> DECIDABLE
  - "A becomes functionally B-like"         -> NOT DECIDABLE (lineage does not measure function)
  - "A and B are the same/different lineage" -> DECIDABLE
  - a claim about a REVERSIBLE STATE within one lineage -> usually NOT decidable, because both
    answers are consistent with the same clonal structure. Think carefully and say why.
  Answer TRUE/FALSE with the reason.

B2 both_answers_actionable  <- THE MOST IMPORTANT FILTER
  Write out, separately and concretely:
    - if the answer is YES (they ARE clonally related): what is the very next experiment?
    - if the answer is NO (they are NOT clonally related): what is the very next experiment?
  If those two are the same thing, or if one of them is vacuous ("we would conclude X and move
  on"), then both_answers_actionable = FALSE and the claim is OUT. Be strict. Most claims fail
  here and that is the expected outcome, not a problem to be avoided.

B3 compartments_obtainable
  Lineage adjudication needs compartment A and compartment B SIMULTANEOUSLY, ideally from the
  same individual. Enumerate exactly which compartments are needed, and state whether they can
  be obtained paired from one human (or one animal). Note explicitly if one compartment is
  only obtainable post-mortem, only at surgery, only during pregnancy, or only in mouse.

B4 already_traced  <- run REAL searches, do not answer from memory
  You MUST run searches covering at least:
    mtDNA + (lineage OR clonal) + <this claim's cell types>
    mtscATAC ; MAESTER ; "single-cell mitochondrial" + <cell type>
    barcode + NK cell ; lentiviral barcoding + lymphocyte
    and for mouse claims: fate mapping ; confetti ; lineage tracing ; Polylox ; CRISPR recorder
  Record TRUE (someone has done it - give PMID and what exactly they did),
  FALSE (searched thoroughly, nobody has), or NOT_SEARCHED.
  RULE G3: "I could not find it" is NOT "it does not exist". If your searches were thin or
  ambiguous, record NOT_SEARCHED rather than FALSE. Being honest here matters more than being
  decisive: a FALSE that should have been NOT_SEARCHED corrupts the shortlist.

RULE G4: if the claim invokes a general process (differentiation, plasticity, residency,
exhaustion, memory), also search whether the T cell / ILC / myeloid fields have already
answered the SAME structural question with lineage tools. A borrowed mechanism's occupancy
often lives in the field it was borrowed from.

TOOLS: python3 /home/user/aaa/scripts/lit.py search "<query>" --n 15
       python3 /home/user/aaa/scripts/lit.py methods <PMID|PMCID> --tag <CLAIM_ID>
       plus WebSearch / WebFetch.`

const B_SCHEMA = {
  type: 'object',
  properties: {
    claim_id: { type: 'string' },
    lineage_decidable: { type: 'string', enum: ['TRUE', 'FALSE'] },
    why: { type: 'string' },
    action_if_yes: { type: 'string' },
    action_if_no: { type: 'string' },
    both_answers_actionable: { type: 'string', enum: ['TRUE', 'FALSE'] },
    both_answers_reason: { type: 'string' },
    compartments_needed: { type: 'string' },
    compartments_obtainable: { type: 'string' },
    compartment_notes: { type: 'string' },
    mouse_only: { type: 'string', enum: ['TRUE', 'FALSE'] },
  },
  required: ['claim_id', 'lineage_decidable', 'why', 'action_if_yes', 'action_if_no',
    'both_answers_actionable', 'both_answers_reason', 'compartments_needed',
    'compartments_obtainable', 'mouse_only'],
}

const OCC_SCHEMA = {
  type: 'object',
  properties: {
    claim_id: { type: 'string' },
    already_traced: { type: 'string', enum: ['TRUE', 'FALSE', 'NOT_SEARCHED'] },
    already_traced_evidence: { type: 'string' },
    nearest_occupant: { type: 'string', description: 'Closest published work with PMID, and precisely how it differs from this claim' },
    crossfield_occupancy: { type: 'string', description: 'Whether the T cell / ILC / myeloid fields have answered the same structural question with lineage tools' },
    search_terms_used: { type: 'array', items: { type: 'string' } },
    search_confidence: { type: 'string', enum: ['HIGH', 'MEDIUM', 'LOW'] },
  },
  required: ['claim_id', 'already_traced', 'already_traced_evidence', 'nearest_occupant',
    'crossfield_occupancy', 'search_terms_used', 'search_confidence'],
}

phase('Decide')

const out = await pipeline(
  CLAIMS,
  (c) => agent(
    `${PREAMBLE}

===== CLAIM =====
ID: ${c.claim_id}
CLAIM: ${c.claim_text}
Phase A strongest evidence: ${c.evidence_type_strongest}
Phase A load-bearing: ${c.load_bearing}
Phase A notes: ${c.notes || 'none'}

Answer B1, B2 and B3. Be strict on B2.`,
    { label: `decide:${c.claim_id}`, phase: 'Decide', schema: B_SCHEMA }
  ),
  (b, c) => b ? agent(
    `${PREAMBLE}

===== B4 OCCUPANCY SEARCH FOR CLAIM ${c.claim_id} =====
CLAIM: ${c.claim_text}
Compartments this claim needs: ${b.compartments_needed}

Run the B4 searches listed above. Your ONLY job is: has anyone already adjudicated this
particular claim with a clonal/lineage method? Report TRUE / FALSE / NOT_SEARCHED, the
nearest occupant with PMID and exactly how it differs, and whether the T cell / ILC / myeloid
fields have already answered the same structural question with lineage tools.

Be honest about search depth. G3: not found != does not exist.`,
    { label: `occupancy:${c.claim_id}`, phase: 'Occupancy', schema: OCC_SCHEMA }
  ).then(o => ({ b, o, claim_id: c.claim_id })) : null
)

const rows = out.filter(Boolean)
log(`Phase B: ${rows.length}/${CLAIMS.length} claims resolved`)
return { n: rows.length, rows }
