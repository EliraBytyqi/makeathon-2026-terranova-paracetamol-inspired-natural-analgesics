import path from 'node:path';
import { pathToFileURL } from 'node:url';

const ROOT = 'C:/Users/HP/makeathon-2026-terranova-paracetamol-inspired-natural-analgesics';
const SKILL = 'C:/Users/HP/.codex/plugins/cache/openai-primary-runtime/presentations/26.1105.12137/skills/presentations';
const { finalizePresentation } = await import(pathToFileURL(path.join(SKILL, 'container_tools/artifact_tool_utils.mjs')).href);
const stagingDir = path.join(ROOT, '.terranova-deck-build');
const result = await finalizePresentation({
  explicitTotalSlideCount: 10,
  requiredNativeTableOwnerSlides: [],
  requiredNativeChartOwnerSlides: [],
  workspaceDir: ROOT,
  candidatePath: path.join(stagingDir, 'candidate.pptx'),
  finalPath: path.join(ROOT, 'deliverables', 'Terranova_Pipeline_Explainer.pptx'),
  pythonExecutable: 'C:/Users/HP/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',
  integrityValidatorPath: path.join(SKILL, 'container_tools/inspect_presentation_package_integrity.py'),
  layoutValidatorPath: path.join(SKILL, 'container_tools/inspect_presentation_layout_geometry.py'),
  layoutArgs: ['--expected-slide-size-emu', '12192000,6858000', '--validate-bullet-geometry', '--validate-heading-fit'],
  fontPolicy: { basis: 'design', families: ['Arial'] },
  verifyArtifactToolImport: true,
  receiptPath: path.join(stagingDir, 'Terranova_Pipeline_Explainer.validation.json'),
});
console.log(JSON.stringify(result));
