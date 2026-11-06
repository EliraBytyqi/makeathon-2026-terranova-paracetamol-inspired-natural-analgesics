/* Domain terms used in this project, each mapped to its ML equivalent.
 *
 * Audience: an ML engineer with no chemistry or biology background. The rule
 * for every entry is that the "equivalent" column must name something the
 * reader already uses, not just restate the definition in simpler words.
 */
window.GLOSSARY = [
  /* ------------------------------------------------- molecules & features */
  {
    term: "SMILES",
    cat: "representation",
    plain: "A molecule written as a text string, e.g. <code>CC(=O)Nc1ccc(O)cc1</code> is paracetamol.",
    ml: "A serialization format for a graph — like writing a DAG as a string. Our primary key for passing molecules between stages.",
  },
  {
    term: "InChIKey",
    cat: "representation",
    plain: "A fixed 27-character identifier derived from the structure itself.",
    ml: "A content hash. Used for deduplication and for the join between datasets, because the same molecule has many valid SMILES spellings.",
  },
  {
    term: "ECFP4 / Morgan fingerprint",
    cat: "representation",
    plain: "A 2048-bit vector where each bit means 'this molecule contains this local substructure'.",
    ml: "Feature hashing over subgraphs — a bag-of-n-grams for molecules. Radius 2 means it looks two bonds out from each atom. This is the model input.",
  },
  {
    term: "Descriptors (MW, cLogP, TPSA, QED)",
    cat: "representation",
    plain: "Scalar physical properties: mass, greasiness, polar surface area, a drug-likeness score.",
    ml: "Hand-engineered numeric features, concatenated onto the fingerprint.",
  },
  {
    term: "Bemis-Murcko scaffold",
    cat: "representation",
    plain: "The molecule's core skeleton with all the decorative side chains stripped off.",
    ml: "<strong>A group key.</strong> This is what the grouped train/test split is grouped on — the exact analogue of splitting by <code>user_id</code> so the same user cannot appear on both sides.",
  },
  {
    term: "Tanimoto similarity",
    cat: "representation",
    plain: "How much two molecules' fingerprints overlap, from 0 to 1.",
    ml: "Jaccard index on binary vectors. Computed as a matrix product, which is why the screen runs in minutes rather than hours.",
  },

  /* ----------------------------------------------------- pharmacology */
  {
    term: "Target",
    cat: "pharmacology",
    plain: "The protein you want the drug to act on. Here: TRPV1, CB1 and FAAH.",
    ml: "One binary classification problem per target. Three targets = three independent label columns, not a multiclass problem — a molecule can hit all three.",
  },
  {
    term: "Ligand",
    cat: "pharmacology",
    plain: "A molecule that binds to a target protein.",
    ml: "A positive example.",
  },
  {
    term: "Agonist / antagonist / inhibitor",
    cat: "pharmacology",
    plain: "Molecules that switch a protein on / block it / shut it down.",
    ml: "The <em>direction</em> of the effect. <strong>We do not model this</strong> — our label is only 'acts on it or not'. A real limitation: an agonist and an antagonist both count as active here.",
  },
  {
    term: "IC50, EC50, Ki, Kd",
    cat: "pharmacology",
    plain: "Four ways of measuring how much of a compound you need for an effect. Lower = more potent.",
    ml: "The continuous target variable, before thresholding. They are not strictly interchangeable, which is a known source of label noise.",
  },
  {
    term: "nM / µM (nanomolar / micromolar)",
    cat: "pharmacology",
    plain: "Concentration units. 1 µM = 1000 nM. <strong>Lower means more potent</strong>, which is the opposite of the intuition most people start with.",
    ml: "Units of the raw label. Everything is normalised to nM before the log transform.",
  },
  {
    term: "pActivity",
    cat: "pharmacology",
    plain: "A rescaling: <code>9 − log10(nM)</code>. Higher = more potent.",
    ml: "A log transform that makes the target roughly linear and flips the sign so bigger is better. 1 µM becomes exactly 6.0.",
  },
  {
    term: "Assay",
    cat: "pharmacology",
    plain: "The lab experiment that produces one measurement.",
    ml: "The labelling process. Different assay types are effectively different annotators with different biases.",
  },
  {
    term: "Binding (B) vs functional (F) assay",
    cat: "pharmacology",
    plain: "Does the molecule stick to the protein, versus does it actually change what the protein does.",
    ml: "Two different label-generation processes measuring related but distinct things. We keep both — dropping functional assays would discard most TRPV1 data.",
  },
  {
    term: "Organism (human / rat / mouse)",
    cat: "pharmacology",
    plain: "Which species the protein was taken from. The same protein differs slightly between species.",
    ml: "Covariate shift between data sources. Kept as an explicit column and modelled as two scopes rather than pooled silently.",
  },
  {
    term: "Natural product",
    cat: "pharmacology",
    plain: "A molecule a living thing makes — here, specifically plants.",
    ml: "The inference set. Structurally very different from the training set, which is the central difficulty of the whole project.",
  },

  /* ------------------------------------------------------------- data */
  {
    term: "ChEMBL",
    cat: "data",
    plain: "A large public database of measured compound–protein activities, curated from published papers.",
    ml: "The labelled training set. 58,394 raw rows in, 31,271 out after curation.",
  },
  {
    term: "LOTUS",
    cat: "data",
    plain: "A database of natural products, each linked to the organism it was found in and a literature DOI.",
    ml: "The unlabelled inference set, 145k molecules. The DOI per row is free provenance — rare and worth showing off.",
  },
  {
    term: "TDC DILI benchmark",
    cat: "data",
    plain: "475 marketed drugs labelled for whether they cause liver injury.",
    ml: "A second, tiny labelled set for the safety model. 462 usable rows is small enough that single-split evaluation is meaningless, hence 5 repeats.",
  },
  {
    term: "Decoy",
    cat: "data",
    plain: "A molecule assumed inactive because nothing suggests it is active.",
    ml: "<strong>Negative sampling.</strong> Real hit rates are well under 1%, so a random drug-like molecule is almost certainly a negative. Labelled <code>presumed_inactive</code>, never <code>inactive</code>, so the assumption stays visible.",
  },
  {
    term: "Property matching (DUD-E style)",
    cat: "data",
    plain: "Choosing decoys with similar mass and greasiness to the real actives.",
    ml: "Negative sampling with covariate balancing. Without it the model learns the shortcut 'heavy and greasy = active' and the metric inflates while transferring nothing.",
  },
  {
    term: "Censored value",
    cat: "data",
    plain: "A measurement that hit the limit of the experiment, recorded as '&gt;10000 nM'.",
    ml: "Right-censored data. Usable as a negative, not as a numeric target. <code>&lt;</code> records are dropped entirely because they cannot confirm a positive.",
  },
  {
    term: "Publication bias",
    cat: "data",
    plain: "Papers report compounds that worked; failures go unpublished.",
    ml: "<strong>The biggest data problem here.</strong> Raw positive rate is 84% for TRPV1, which is not a base rate. This is why decoys exist.",
  },
  {
    term: "Ambiguous band (1–10 µM)",
    cat: "data",
    plain: "Potency between the two thresholds — too weak to call active, too strong to call inactive.",
    ml: "Deliberately dropped rows. Forcing a label on genuinely borderline examples injects noise; roughly 1,800 CB1 rows are discarded this way.",
  },
  {
    term: "Leakage control",
    cat: "data",
    plain: "Removing library compounds that already appear in the training data.",
    ml: "Standard train/inference deduplication, on InChIKey. 216 compounds removed and kept aside as a sanity check.",
  },

  /* ------------------------------------------------------- evaluation */
  {
    term: "Scaffold split",
    cat: "evaluation",
    plain: "Train and test share no core molecular skeleton.",
    ml: "<code>GroupShuffleSplit</code> on scaffold. An out-of-distribution evaluation. <strong>This is the honest number</strong>, because plant chemistry shares almost no scaffolds with the training set.",
  },
  {
    term: "Random split",
    cat: "evaluation",
    plain: "An ordinary shuffle.",
    ml: "The optimistic baseline. It leaks: chemical series are families of near-identical molecules, so a shuffle puts siblings on both sides and the model memorises the family. Reported only for contrast.",
  },
  {
    term: "PR-AUC",
    cat: "evaluation",
    plain: "Area under the precision–recall curve.",
    ml: "The lead metric. Correct choice for imbalanced retrieval — ROC-AUC looks great here regardless, because true negatives dominate.",
  },
  {
    term: "EF 1% (enrichment factor)",
    cat: "evaluation",
    plain: "How many more actives you find in the top 1% than by picking at random.",
    ml: "<strong>Lift@k.</strong> The literal deployment metric: we hand over a top-50 list, so this is the only number that reflects what actually happens.",
  },
  {
    term: "MCC",
    cat: "evaluation",
    plain: "Matthews correlation coefficient.",
    ml: "A balanced single-number classification score, robust to class imbalance. Reported at a fixed 0.5 threshold.",
  },
  {
    term: "Applicability domain (AD)",
    cat: "evaluation",
    plain: "How similar a molecule is to the closest thing the model was trained on.",
    ml: "<strong>OOD / novelty detection.</strong> Max Tanimoto to the training set. Below 0.25 the prediction is extrapolation. 20.8% of the library falls outside it.",
  },
  {
    term: "Composite score",
    cat: "evaluation",
    plain: "The final ranking number combining activity, safety, novelty and brain penetration.",
    ml: "A hand-weighted reranker over the classifier output. Weights are hardcoded and printed rather than tuned, because there is no held-out set to tune them against.",
  },

  /* ----------------------------------------------------------- safety */
  {
    term: "Metabolite",
    cat: "safety",
    plain: "What your body chemically converts a drug into after you take it.",
    ml: "A hidden transformation between the input you model and the thing that actually acts. The molecule you score is not always the molecule that matters.",
  },
  {
    term: "Reactive metabolite",
    cat: "safety",
    plain: "A metabolite that is chemically aggressive and attacks proteins in the cell.",
    ml: "The hazard being filtered. Not visible in any activity label — it needs its own rule layer.",
  },
  {
    term: "NAPQI",
    cat: "safety",
    plain: "The specific reactive metabolite that destroys livers in paracetamol overdose.",
    ml: "The concrete thing this project filters against. Severity-3 alerts are patterns that can form it.",
  },
  {
    term: "SMARTS",
    cat: "safety",
    plain: "A pattern language for matching substructures inside molecules.",
    ml: "<strong>Regex, for molecular graphs.</strong> The 25 structural alerts are 25 patterns.",
  },
  {
    term: "Structural alert",
    cat: "safety",
    plain: "A substructure associated with a known hazard.",
    ml: "A rule-based filter running alongside the learned model. It catches what the training data cannot, which is the entire argument for keeping it.",
  },
  {
    term: "Quinone imine / Michael acceptor / catechol / epoxide",
    cat: "safety",
    plain: "Families of chemical groups that are prone to attacking proteins.",
    ml: "The specific patterns the alert rules match. You do not need to know the chemistry — only that each has documented precedent, and a severity from 1 to 3.",
  },
  {
    term: "PAINS",
    cat: "safety",
    plain: "Substructures notorious for producing fake positives in lab assays.",
    ml: "<strong>Known label-noise patterns.</strong> Molecules containing them test positive against everything. 13,276 library compounds are flagged and excluded from the shortlist.",
  },
  {
    term: "DILI",
    cat: "safety",
    plain: "Drug-induced liver injury.",
    ml: "The binary label of the fourth classifier. Notoriously hard from structure alone, because much clinical risk is dose- and patient-dependent.",
  },
  {
    term: "CNS-like",
    cat: "safety",
    plain: "Likely to cross from blood into the brain.",
    ml: "A heuristic boolean feature from mass/greasiness/polarity thresholds. Relevant because the hypothesis is <em>central</em> pain relief, so a compound that cannot reach the brain is a weaker candidate.",
  },
  {
    term: "Severity 1 / 2 / 3",
    cat: "safety",
    plain: "How strong the evidence is that an alert leads to a real hazard.",
    ml: "Rule priority. 3 = direct documented precedent for the paracetamol mechanism; 1 = context-dependent. Enters the ranking as an exponential penalty.",
  },
];
