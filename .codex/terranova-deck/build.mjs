import fs from 'node:fs/promises';
import path from 'node:path';
import { Presentation, PresentationFile } from 'file:///C:/Users/HP/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';

const ROOT = 'C:/Users/HP/makeathon-2026-terranova-paracetamol-inspired-natural-analgesics';
const OUT = path.join(ROOT, 'deliverables');
const PREVIEW = path.join(ROOT, '.codex/terranova-deck/previews');
const C = {
  ink: '#173C35', pine: '#244F45', sage: '#AFC9B4', mint: '#DCE9DC',
  paper: '#F7F4EA', white: '#FFFFFF', muted: '#65756E', line: '#D6DED4',
  orange: '#D9794B', paleOrange: '#F5E2D6', blue: '#4C8292', paleBlue: '#DFEBEF',
  gold: '#D4AD5F', dark: '#102B27', red: '#A85142', paleRed: '#F2DDDA',
};
const W = 1280, H = 720;
const ppt = Presentation.create({ slideSize: { width: W, height: H } });

function rect(s, x, y, w, h, fill, radius = false, line = 'none') {
  return s.shapes.add({ geometry: radius ? 'roundRect' : 'rect', position: { left: x, top: y, width: w, height: h }, fill, line: { style: 'solid', fill: line, width: line === 'none' ? 0 : 1 }, ...(radius ? { borderRadius: 'rounded-xl' } : {}) });
}
function text(s, str, x, y, w, h, size = 24, color = C.ink, bold = false, align = 'left') {
  const z = s.shapes.add({ geometry: 'textbox', position: { left: x, top: y, width: w, height: h }, fill: 'none', line: { style: 'solid', fill: 'none', width: 0 } });
  z.text = str;
  z.text.style = { fontFamily: 'Aptos', fontSize: size, color, bold, alignment: align, verticalAlignment: 'middle' };
  return z;
}
function rule(s, x, y, w, color = C.line, h = 2) { rect(s, x, y, w, h, color); }
function base(title, n, dark = false) {
  const s = ppt.slides.add();
  s.background.fill = dark ? C.dark : C.paper;
  if (title) text(s, title, 68, 40, 1120, 58, 40, dark ? C.paper : C.ink, true);
  rule(s, 68, 112, 1144, dark ? '#36534B' : C.line, 2);
  text(s, 'TERRANOVA  /  MAKEATHON 2026', 68, 675, 700, 20, 13, dark ? '#AEC4B7' : C.muted, true);
  text(s, String(n).padStart(2, '0'), 1150, 672, 62, 24, 15, dark ? '#AEC4B7' : C.muted, true, 'right');
  return s;
}
function note(s, refs, extra = '') { s.speakerNotes.text = `${extra}${extra ? '\n\n' : ''}Sources: ${refs}`; }
function label(s, value, x, y, w, color = C.pine) { text(s, value.toUpperCase(), x, y, w, 24, 14, color, true); }
function flowNode(s, x, y, w, num, heading, detail, color = C.pine) {
  text(s, num, x, y, 48, 38, 24, color, true);
  text(s, heading, x + 52, y, w - 52, 36, 22, C.ink, true);
  text(s, detail, x + 52, y + 40, w - 58, 72, 18, C.muted, false);
}

// 1 — Cover
{
  const s = ppt.slides.add(); s.background.fill = C.dark;
  // Quiet botanical-inspired composition built from editable shapes.
  rect(s, 918, 0, 362, 720, '#183B34');
  rect(s, 1028, 70, 5, 565, '#557C64');
  for (const [x,y,w,h,c] of [[864,130,220,64,'#355F4B'],[1000,240,210,62,'#466E52'],[836,350,205,64,'#2B5947'],[1004,468,196,60,'#466E52'],[868,564,190,56,'#355F4B']]) {
    const sh = rect(s,x,y,w,h,c,true); sh.rotation = (x < 950 ? -28 : 28);
  }
  text(s, 'TERRANOVA', 76, 84, 740, 82, 58, C.paper, true);
  text(s, 'Plant compounds for the\nparacetamol analgesic pathway', 76, 194, 730, 154, 37, '#DCE9DC', true);
  rule(s, 78, 382, 166, C.orange, 6);
  text(s, 'A computational discovery pipeline', 78, 414, 700, 44, 25, C.paper);
  text(s, 'Makeathon 2026 · Track 02 · Natural Products', 78, 494, 700, 32, 18, '#B6C7BB');
  text(s, 'Elira Bytyqi   /   Shkurta Jashari   /   Fjolla Nebihi', 78, 540, 770, 34, 18, '#B6C7BB');
  text(s, 'COMPUTATIONAL PREDICTIONS ONLY', 78, 640, 640, 24, 14, '#E4B18C', true);
  note(s, 'Project README.md; results/REPORT.md', 'Terranova screens hypotheses for laboratory follow-up. No compound was assayed by the team.');
}

// 2 — Biological rationale
{
  const s = base('The biological idea', 2);
  text(s, 'Paracetamol metabolism points to two chemically distinct routes.', 70, 132, 1120, 36, 23, C.muted);
  rect(s, 72, 198, 534, 314, C.paleBlue, true);
  rect(s, 672, 198, 534, 314, C.paleOrange, true);
  label(s, 'Pain pathway', 102, 222, 220, C.blue);
  text(s, 'Paracetamol', 102, 264, 180, 40, 24, C.ink, true);
  text(s, 'p-aminophenol', 290, 264, 220, 40, 22, C.ink, true);
  text(s, 'FAAH', 178, 330, 100, 34, 19, C.blue, true, 'center');
  text(s, 'AM404', 330, 330, 130, 34, 24, C.ink, true, 'center');
  rule(s, 252, 299, 54, C.blue, 3); rule(s, 272, 295, 12, C.blue, 11);
  rule(s, 296, 346, 24, C.blue, 3); rule(s, 310, 342, 12, C.blue, 11);
  text(s, 'Associated with TRPV1 and cannabinoid signalling involving CB1', 102, 402, 424, 68, 20, C.ink);
  label(s, 'Overdose toxicity route', 702, 222, 310, C.red);
  text(s, 'Paracetamol', 702, 264, 180, 40, 24, C.ink, true);
  text(s, 'CYP enzymes', 895, 264, 200, 40, 22, C.ink, true);
  text(s, 'NAPQI', 805, 330, 150, 34, 25, C.red, true, 'center');
  text(s, 'Glutathione depletion can allow liver-cell damage', 702, 402, 444, 68, 20, C.ink);
  rule(s, 861, 299, 24, C.red, 3); rule(s, 875, 295, 12, C.red, 11);
  rect(s, 72, 547, 1134, 74, C.ink, true);
  text(s, 'Research question', 96, 560, 200, 28, 16, '#BFD3C2', true);
  text(s, 'Can plant chemistry offer testable hypotheses around the analgesic targets while flagging potential liabilities?', 96, 589, 1070, 30, 21, C.white, true);
  note(s, 'Project README.md; results/REPORT.md; pitch script supplied by the team', 'The pathways motivate the screen. A predicted target interaction does not prove analgesia or mechanism.');
}

// 3 — End-to-end pipeline
{
  const s = base('Terranova pipeline', 3);
  text(s, 'A traceable path from public evidence to experiments worth considering.', 70, 132, 1120, 36, 23, C.muted);
  const xs = [76, 314, 552, 790, 1028];
  const heads = ['Collect', 'Curate', 'Train', 'Screen', 'Prioritise'];
  const det = ['ChEMBL activity\nLOTUS plants', 'Standardise\nstructures + labels', 'TRPV1 · CB1 · FAAH\nRF + LightGBM', '145,247 plant\nstructures scored', 'Safety + similarity\n50 hypotheses'];
  const cols = [C.blue, C.sage, C.gold, C.orange, C.pine];
  for (let i=0; i<5; i++) {
    rect(s, xs[i], 216, 170, 12, cols[i]);
    text(s, String(i+1).padStart(2,'0'), xs[i], 249, 48, 35, 22, cols[i], true);
    text(s, heads[i], xs[i], 293, 205, 39, 24, C.ink, true);
    text(s, det[i], xs[i], 343, 195, 90, 18, C.muted);
    if (i<4) { rule(s, xs[i]+177, 221, 45, C.line, 2); }
  }
  rule(s, 78, 466, 1124, C.line, 2);
  text(s, 'What users receive', 78, 492, 250, 32, 21, C.ink, true);
  text(s, 'A ranked candidate table with predicted target scores, applicability estimates, safety flags, plant provenance and literature references.', 78, 532, 1050, 62, 22, C.muted);
  note(s, 'README.md; app/README.md; results/REPORT.md');
}

// 4 — Evidence and datasets
{
  const s = base('Evidence behind the screen', 4);
  text(s, 'Three data roles keep learning, discovery and evaluation distinct.', 70, 132, 1120, 34, 22, C.muted);
  const x = [76, 465, 854], widths = 344;
  label(s, 'Training evidence', x[0], 198, 240, C.blue);
  text(s, 'ChEMBL', x[0], 237, widths[0], 40, 29, C.ink, true);
  text(s, '58,394 raw activity records', x[0], 294, widths[0], 34, 21, C.ink, true);
  text(s, '31,271 records after curation\nTRPV1 · CB1 · FAAH', x[0], 338, widths[0], 70, 19, C.muted);
  rect(s, x[0], 424, 282, 8, C.blue);
  label(s, 'Natural-product library', x[1], 198, 300, C.pine);
  text(s, 'LOTUS', x[1], 237, widths[1], 40, 29, C.ink, true);
  text(s, '145,380 plant structures', x[1], 294, widths[1], 34, 21, C.ink, true);
  text(s, '133 exact overlaps removed\n145,247 screened', x[1], 338, widths[1], 70, 19, C.muted);
  rect(s, x[1], 424, 282, 8, C.sage);
  label(s, 'External evaluation', x[2], 198, 260, C.orange);
  text(s, 'PubChem assays', x[2], 237, widths[2], 40, 29, C.ink, true);
  text(s, 'Measured outcomes', x[2], 294, widths[2], 34, 21, C.ink, true);
  text(s, 'Expanded holdout set\nnot used to train models', x[2], 338, widths[2], 70, 19, C.muted);
  rect(s, x[2], 424, 282, 8, C.orange);
  rect(s, 76, 486, 1128, 115, C.mint, true);
  text(s, 'Curation matters', 102, 505, 250, 32, 22, C.ink, true);
  text(s, 'Structures and measurements were standardised. Ambiguous potency records were excluded from class labels. Exact training overlaps were removed from the screen to reduce leakage.', 102, 545, 1045, 45, 19, C.ink);
  note(s, 'results/REPORT.md; results/pubchem_validation_expanded.json; results/screening_stats.json');
}

// 5 — Training data challenge
{
  const s = base('Solving the label imbalance problem', 5);
  text(s, 'Published activity data favour compounds that researchers already found interesting.', 70, 132, 1120, 34, 22, C.muted);
  text(s, '84%', 82, 208, 250, 120, 76, C.orange, true);
  text(s, 'of curated TRPV1 measurements\nwere labelled active', 88, 328, 395, 68, 23, C.ink, true);
  rule(s, 500, 202, 2, C.line, 352);
  flowNode(s, 550, 210, 590, '01', 'Matched presumed-inactive decoys', 'Balance training examples while keeping assumed negatives distinct from measured inactive compounds.', C.orange);
  flowNode(s, 550, 328, 590, '02', 'Lower weight for decoys', 'Presumed inactives receive lower training weight, so the model treats them as less certain.', C.blue);
  flowNode(s, 550, 446, 590, '03', 'Scaffold-based evaluation', 'Hold out chemical scaffolds to test performance on structures beyond close analogues.', C.pine);
  rect(s, 82, 506, 397, 78, C.paleOrange, true);
  text(s, 'Decoys are presumed inactive, not experimentally confirmed negatives.', 102, 518, 350, 54, 18, C.red, true);
  note(s, 'results/REPORT.md; data/processed/dataset_summary.csv; results/training_report.json', 'The 84% figure refers to curated TRPV1 measurements. Decoy sample weight is 0.4 per the project report.');
}

// 6 — The two activity model families
{
  const s = base('Two model families score target activity', 6);
  text(s, 'Each target model sees the same molecular representation and predicts an activity ranking.', 70, 132, 1120, 36, 22, C.muted);
  // Editable molecule feature representation
  rect(s, 78, 202, 267, 308, C.mint, true);
  label(s, 'Molecular input', 104, 226, 210, C.pine);
  text(s, '2,048-bit Morgan\nfingerprint', 104, 268, 212, 66, 25, C.ink, true);
  text(s, 'Circular substructure\nfeatures', 104, 344, 212, 52, 18, C.muted);
  rule(s, 104, 412, 208, C.line, 2);
  text(s, '+ 14 molecular\nproperty descriptors', 104, 430, 212, 58, 20, C.ink, true);
  rect(s, 410, 202, 340, 308, C.paleBlue, true);
  label(s, 'Random Forest', 438, 226, 260, C.blue);
  text(s, 'Many decision trees\ncombine their votes', 438, 274, 280, 75, 26, C.ink, true);
  text(s, 'Good baseline for nonlinear structure–activity patterns.', 438, 379, 280, 75, 19, C.muted);
  rect(s, 790, 202, 340, 308, C.paleOrange, true);
  label(s, 'LightGBM', 818, 226, 260, C.orange);
  text(s, 'Trees are added\nsequentially', 818, 274, 280, 75, 26, C.ink, true);
  text(s, 'Each round focuses on errors left by earlier trees.', 818, 379, 280, 75, 19, C.muted);
  text(s, 'For TRPV1, CB1 and FAAH, the team compared both using scaffold splits and reported model-specific metrics.', 80, 548, 1090, 54, 21, C.ink, true);
  note(s, 'src/featurize.py; src/train.py; results/model_metrics.csv; results/REPORT.md', 'The two classifiers are Random Forest and LightGBM. Output scores are model predictions, not confirmed probabilities of efficacy. Activity labels do not distinguish agonists from antagonists.');
}

// 7 — Safety and applicability
{
  const s = base('Safety context and applicability domain', 7);
  text(s, 'A high target score is only one part of a candidate’s evidence.', 70, 132, 1120, 34, 22, C.muted);
  rect(s, 78, 208, 510, 344, C.paleOrange, true);
  label(s, 'Chemical alerts', 108, 238, 280, C.red);
  text(s, '25 structural patterns', 108, 280, 395, 40, 29, C.ink, true);
  text(s, 'Flag motifs associated with reactive metabolites, including para-aminophenol and quinone chemistry.', 108, 338, 420, 88, 20, C.ink);
  text(s, 'Mechanism-based flags help catch chemistry a general toxicity model can miss.', 108, 452, 420, 63, 18, C.muted);
  rect(s, 638, 208, 560, 344, C.paleBlue, true);
  label(s, 'Hepatic-liability panel', 668, 238, 340, C.blue);
  text(s, 'Three evidence heads', 668, 280, 440, 40, 29, C.ink, true);
  text(s, 'Clinical DILI labels\nMitochondrial membrane-potential data\nHepG2 cytotoxicity data', 668, 338, 470, 116, 20, C.ink);
  text(s, 'The model is a weak prior. It cannot establish safety and may miss paracetamol’s dose-dependent toxicity.', 668, 462, 470, 58, 18, C.muted);
  rect(s, 78, 578, 1120, 48, C.ink, true);
  text(s, 'Applicability domain = chemical familiarity estimate, not probability of activity.', 100, 586, 1060, 32, 20, C.white, true, 'center');
  note(s, 'src/safety.py; src/train_hepatic.py; results/REPORT.md; results/hepatic_report.json');
}

// 8 — Validation results
{
  const s = base('External evaluation shows uneven performance', 8);
  text(s, 'ROC-AUC measures how well scores rank measured actives above measured inactives. Chance = 0.50.', 70, 132, 1120, 36, 20, C.muted);
  const entries = [
    { name:'TRPV1', val:0.918, color:C.blue, msg:'Strong overall; novel-space check has too few inactives for the strictest slice.' },
    { name:'CB1', val:0.917, color:C.orange, msg:'Strong overall; falls to 0.711 in the most dissimilar evaluated chemistry.' },
    { name:'FAAH', val:0.625, color:C.pine, msg:'Weak ranking performance; improve evidence coverage before relying on it.' },
  ];
  entries.forEach((d,i)=>{
    const y=218+i*112;
    text(s,d.name,82,y,130,38,24,C.ink,true);
    rect(s,230,y+7,510,26,'#E1E4DA',true);
    rect(s,230,y+7,510*d.val,26,d.color,true);
    text(s,d.val.toFixed(3),760,y-2,100,40,24,C.ink,true);
    text(s,d.msg,870,y-2,330,66,17,C.muted);
  });
  rule(s,230,567,510,'#B8C3B8',2);
  text(s,'0.50 chance',220,579,150,23,14,C.muted);
  text(s,'1.00 perfect ranking',610,579,175,23,14,C.muted,false,'right');
  rect(s,82,622,1110,34,C.paleOrange,true);
  text(s,'Separate literature check: 9 held-out compounds · median percentile 55.1 · not statistically significant · inconclusive',95,625,1085,26,16,C.red,true,'center');
  note(s, 'results/pubchem_validation_expanded.json; results/literature_validation_summary.json; pitch script supplied by the team', 'PubChem expanded validation ROC-AUC values: TRPV1 0.9176, CB1 0.9174, FAAH 0.6245. Strict novel-space CB1 ROC-AUC below 0.5 similarity: 0.7108. The literature validation is a separate, small test and is inconclusive.');
}

// 9 — Product outputs
{
  const s = base('From predictions to an inspectable shortlist', 9);
  text(s, 'A candidate record keeps the prediction attached to its context.', 70, 132, 1120, 36, 22, C.muted);
  const rows = [
    ['33','TRPV1','prioritised hypotheses',C.blue],
    ['15','CB1','prioritised hypotheses',C.orange],
    ['2','FAAH','prioritised hypotheses',C.pine],
  ];
  rows.forEach((r,i)=>{
    const y=220+i*108;
    text(s,r[0],86,y,100,60,46,r[3],true);
    text(s,r[1],206,y+4,130,38,24,C.ink,true);
    text(s,r[2],344,y+5,290,34,19,C.muted);
    rule(s,84,y+76,544,C.line,2);
  });
  rect(s, 700, 206, 492, 340, C.mint, true);
  label(s,'Dashboard + knowledge graph',730,232,390,C.pine);
  text(s,'Explore',730,276,160,34,24,C.ink,true);
  text(s,'Predicted target scores and applicability estimate',730,318,410,46,18,C.muted);
  text(s,'Inspect',730,374,160,34,24,C.ink,true);
  text(s,'Safety alerts, plant sources and literature DOIs',730,416,410,46,18,C.muted);
  text(s,'Trace',730,472,160,34,24,C.ink,true);
  text(s,'Measured, literature and predicted edges are labelled separately',730,514,420,46,18,C.muted);
  text(s,'The top 50 are hypotheses for experiments. They are not discovered analgesics.',84,585,1070,40,22,C.ink,true);
  note(s, 'results/candidates_top50.csv; results/knowledge_graph.html; results/kg_stats.json; app/README.md', 'The shortlist distribution is 33 TRPV1, 15 CB1 and 2 FAAH. The graph contains 9,090 nodes and 9,846 edges in the reported build.');
}

// 10 — Lessons and next steps
{
  const s = base('What the pipeline solved, and what comes next', 10, true);
  text(s, 'The result is a clearer experiment queue, with uncertainty attached.', 70, 132, 1110, 38, 23, '#B6C7BB');
  text(s, 'What we solved', 78, 216, 430, 36, 25, '#B7D0BA', true);
  text(s, 'Connected public bioactivity data to a large plant library.\n\nReduced imbalance and scaffold leakage risks.\n\nAdded chemical safety flags and traceable evidence.', 78, 270, 490, 220, 22, C.paper);
  rule(s, 640, 208, 2, '#49675A', 372);
  text(s, 'Next experiments', 704, 216, 430, 36, 25, '#E8B18C', true);
  text(s, 'Expand FAAH bioactivity coverage.\n\nTest candidates in target-specific assays.\n\nFollow with selectivity and cytotoxicity studies.', 704, 270, 470, 220, 22, C.paper);
  rect(s, 78, 548, 1120, 72, '#244A3E', true);
  text(s, '50 testable hypotheses   /   no compounds assayed by the team   /   computational predictions only', 104, 566, 1068, 34, 20, C.white, true, 'center');
  text(s, 'Elira Bytyqi · Shkurta Jashari · Fjolla Nebihi', 80, 636, 830, 28, 18, '#B6C7BB');
  note(s, 'README.md; results/REPORT.md; results/literature_validation_summary.json', 'Close by emphasizing that the shortlist proposes what to test next. It does not show efficacy or safety.');
}

await fs.mkdir(OUT, { recursive: true });
await fs.mkdir(PREVIEW, { recursive: true });
for (let i=0; i<ppt.slides.items.length; i++) {
  const slide = ppt.slides.items[i];
  const blob = await ppt.export({ slide, format: 'png', scale: 1 });
  await fs.writeFile(path.join(PREVIEW, `slide-${String(i+1).padStart(2,'0')}.png`), new Uint8Array(await blob.arrayBuffer()));
}
const staging = path.join(ROOT,'.codex/terranova-deck/candidate.pptx');
const draft = await PresentationFile.exportPptx(ppt);
await draft.save(staging);
console.log(JSON.stringify({ slides: ppt.slides.items.length, staging, preview: PREVIEW }));
