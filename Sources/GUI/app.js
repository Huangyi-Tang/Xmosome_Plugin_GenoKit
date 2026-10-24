// GenoKit GUI Controller — Aligned with PedigreeKit Architecture

let bedBytes = null;
let bimText = null;
let famText = null;
let bedFileName = '';
let bimFileName = '';
let famFileName = '';

// DOM Elements
const inputPanel = document.getElementById('inputPanel');
const inputToggle = document.getElementById('inputToggle');
const inputSummary = document.getElementById('inputSummary');

const chooseFileButton = document.getElementById('chooseFileButton');
const clearButton = document.getElementById('clearButton');
const fileInput = document.getElementById('fileInput');
const fileName = document.getElementById('fileName');

const badgeBed = document.getElementById('badgeBed');
const badgeBim = document.getElementById('badgeBim');
const badgeFam = document.getElementById('badgeFam');

const mindEnabledInput = document.getElementById('mindEnabledInput');
const mindInput = document.getElementById('mindInput');
const genoEnabledInput = document.getElementById('genoEnabledInput');
const genoInput = document.getElementById('genoInput');
const hweEnabledInput = document.getElementById('hweEnabledInput');
const hweInput = document.getElementById('hweInput');
const mafEnabledInput = document.getElementById('mafEnabledInput');
const mafInput = document.getElementById('mafInput');
const maxMafEnabledInput = document.getElementById('maxMafEnabledInput');
const maxMafInput = document.getElementById('maxMafInput');

const orderSelect = document.getElementById('orderSelect');
const roundsInput = document.getElementById('roundsInput');
const makeBedInput = document.getElementById('makeBedInput');

const runBtn = document.getElementById('runBtn');
const computeNote = document.getElementById('computeNote');

const resultsPanel = document.getElementById('resultsPanel');
const statSamples = document.getElementById('statSamples');
const statMind = document.getElementById('statMind');
const statVariants = document.getElementById('statVariants');
const statGeno = document.getElementById('statGeno');
const statHwe = document.getElementById('statHwe');
const statMaf = document.getElementById('statMaf');
const downloadButtons = document.getElementById('downloadButtons');
const reportContainer = document.getElementById('reportContainer');

const progressModal = document.getElementById('progressModal');
const progressFill = document.getElementById('progressFill');
const progressPercent = document.getElementById('progressPercent');
const progressStage = document.getElementById('progressStage');
const toastHost = document.getElementById('toastHost');

// Declarative GUI → host synchronization. The host combines these values with
// manifest `cli.arguments` and native file URLs to render the universal preview.
function syncCLIState() {
  const binaryNames = [bedFileName, bimFileName, famFileName].filter(Boolean);
  const binaryPrefixes = binaryNames.map(name => name.replace(/\.(bed|bim|fam)$/i, ''));
  const bfile = binaryPrefixes.length === 3 && binaryPrefixes.every(prefix => prefix === binaryPrefixes[0])
    ? binaryPrefixes[0]
    : null;
  window.xmosomeCLI?.setState({
    arguments: {
      bfile,
      bed: bfile ? null : (bedFileName || null),
      bim: bfile ? null : (bimFileName || null),
      fam: bfile ? null : (famFileName || null),
      mind: optionalThresholdValue(mindEnabledInput, mindInput, 0.1),
      geno: optionalThresholdValue(genoEnabledInput, genoInput, 0.05),
      hwe: optionalThresholdValue(hweEnabledInput, hweInput, 0.000001),
      maf: optionalThresholdValue(mafEnabledInput, mafInput, 0.05),
      max_maf: optionalThresholdValue(maxMafEnabledInput, maxMafInput, 0.4),
      order: orderSelect?.value || 'mind,geno,hwe,maf',
      rounds: Number.parseInt(roundsInput?.value || '1', 10) || 1,
      make_bed: makeBedInput?.checked ?? false,
    }
  });
}

/** Spinner step matching the default value's order of magnitude (0.1 → 0.1, 0.05 → 0.01, 1e-6 → 1e-6). */
function stepForDefault(value) {
  const n = Number(value);
  if (!Number.isFinite(n) || n === 0) return 0.01;
  return Math.pow(10, Math.floor(Math.log10(Math.abs(n))));
}

function optionalThresholdValue(enabledInput, valueInput, fallback) {
  if (!enabledInput?.checked) return null;
  const parsed = Number.parseFloat(valueInput?.value || String(fallback));
  return Number.isFinite(parsed) ? parsed : fallback;
}

function bindOptionalThreshold(enabledInput, valueInput, defaultValue) {
  if (!enabledInput || !valueInput) return;
  const step = stepForDefault(defaultValue);
  valueInput.step = String(step);
  if (!valueInput.value) valueInput.value = String(defaultValue);
  const syncEnabled = () => {
    valueInput.disabled = !enabledInput.checked;
    syncCLIState();
  };
  enabledInput.addEventListener('change', syncEnabled);
  valueInput.addEventListener('input', syncCLIState);
  valueInput.addEventListener('change', syncCLIState);
  syncEnabled();
}

function notice(message, isError = false) {
  showToast(message, isError ? 'error' : 'ok');
}

function base64ToBytes(b64) {
  const raw = atob(b64);
  const arr = new Uint8Array(raw.length);
  for (let i = 0; i < raw.length; i++) arr[i] = raw.charCodeAt(i);
  return arr;
}

function applyNamedContent(name, content, encoding) {
  const lower = (name || '').toLowerCase();
  if (lower.endsWith('.bed')) {
    bedBytes = encoding === 'base64' ? base64ToBytes(content || '') : new TextEncoder().encode(content || '');
    bedFileName = name;
    badgeBed?.classList.add('loaded');
    return true;
  }
  if (lower.endsWith('.bim')) {
    bimText = content || '';
    bimFileName = name;
    badgeBim?.classList.add('loaded');
    return true;
  }
  if (lower.endsWith('.fam')) {
    famText = content || '';
    famFileName = name;
    badgeFam?.classList.add('loaded');
    return true;
  }
  return false;
}

// Host (Mac / iPad) delivers native multi-file picks through this callback.
function handleReceiveFiles(payload) {
  if (!payload?.ok) {
    notice(payload?.error || 'File selection cancelled.', true);
    return;
  }
  const files = Array.isArray(payload.files) ? payload.files : [];
  if (!files.length) {
    notice('No genotype files were selected.', true);
    return;
  }
  let loaded = 0;
  let sawPedMap = false;
  for (const file of files) {
    const name = (file.filename || '').toLowerCase();
    if (name.endsWith('.ped') || name.endsWith('.map')) sawPedMap = true;
    if (applyNamedContent(file.filename || '', file.content || '', file.encoding || 'utf8')) {
      loaded += 1;
    }
  }
  checkReady();
  if (loaded === 0) {
    notice(
      sawPedMap
        ? 'GenoKit GUI loads .bed/.bim/.fam. For .ped/.map use CLI: genokit --file <prefix>'
        : 'Choose .bed / .bim / .fam files (same prefix preferred).',
      true
    );
  } else if (!(bedBytes && bimText && famText)) {
    notice(`Loaded ${loaded} file(s). Still waiting for the full BED/BIM/FAM set.`, false);
  }
}

// Prefer the host stub's impl slot so early deliveries queue correctly.
window.__xmosomeReceiveFilesImpl = handleReceiveFiles;
window.xmosomeReceiveFiles = function (payload) {
  if (typeof window.__xmosomeReceiveFilesImpl === 'function') {
    return window.__xmosomeReceiveFilesImpl(payload);
  }
  window.__xmosomePendingReceiveFiles = payload;
};

// Drain a payload that arrived before this script finished loading.
if (window.__xmosomePendingReceiveFiles) {
  const pending = window.__xmosomePendingReceiveFiles;
  window.__xmosomePendingReceiveFiles = null;
  handleReceiveFiles(pending);
}

// Single-file fallback used by the injected host bridge / older hosts.
window.xmosomePickDataResult = (ok, filename, content, errorMessage) => {
  if (!ok) {
    notice(errorMessage || 'File selection cancelled.', true);
    return;
  }
  if (!applyNamedContent(filename || '', content || '', 'utf8')) {
    notice(errorMessage || 'GenoKit expects .bed/.bim/.fam (native multi-file pick).', true);
    return;
  }
  checkReady();
};

function pickFile() {
  if (window.webkit?.messageHandlers?.xmosomePickData) {
    window.webkit.messageHandlers.xmosomePickData.postMessage({ mode: 'genotype' });
    notice('Choose genotype files…');
    return;
  }
  fileInput?.click();
}

// 1. Collapsible Input Panel Toggle
if (inputToggle && inputPanel) {
  inputToggle.addEventListener('click', e => {
    // Prevent collapsing when clicking on help tip inside header
    if (e.target.closest('.helpTip')) return;
    const isExpanded = inputPanel.classList.toggle('is-expanded');
    inputToggle.setAttribute('aria-expanded', isExpanded ? 'true' : 'false');
  });
  inputToggle.addEventListener('keydown', e => {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      const isExpanded = inputPanel.classList.toggle('is-expanded');
      inputToggle.setAttribute('aria-expanded', isExpanded ? 'true' : 'false');
    }
  });
}

// 2. Help Tip Mobile / Click Accessibility
document.querySelectorAll('.helpTip').forEach(tip => {
  tip.addEventListener('click', e => {
    e.stopPropagation();
    tip.classList.toggle('is-open');
  });
});
document.addEventListener('click', () => {
  document.querySelectorAll('.helpTip.is-open').forEach(t => t.classList.remove('is-open'));
});

// 3. Optional threshold controls (checkbox + number, HWE-style)
bindOptionalThreshold(mindEnabledInput, mindInput, 0.1);
bindOptionalThreshold(genoEnabledInput, genoInput, 0.05);
bindOptionalThreshold(hweEnabledInput, hweInput, 0.000001);
bindOptionalThreshold(mafEnabledInput, mafInput, 0.05);
bindOptionalThreshold(maxMafEnabledInput, maxMafInput, 0.4);

// 4. File Drop & Pick Handling
if (chooseFileButton && fileInput) {
  chooseFileButton.addEventListener('click', pickFile);
  chooseFileButton.addEventListener('dragover', e => {
    e.preventDefault();
    chooseFileButton.classList.add('is-drop-target');
  });
  chooseFileButton.addEventListener('dragleave', () => chooseFileButton.classList.remove('is-drop-target'));
  chooseFileButton.addEventListener('drop', e => {
    e.preventDefault();
    chooseFileButton.classList.remove('is-drop-target');
    if (e.dataTransfer.files?.length) {
      handleFiles(Array.from(e.dataTransfer.files));
    }
  });
  fileInput.addEventListener('change', () => {
    if (fileInput.files?.length) {
      handleFiles(Array.from(fileInput.files));
    }
  });
}

if (clearButton) {
  clearButton.addEventListener('click', () => {
    bedBytes = null;
    bimText = null;
    famText = null;
    bedFileName = '';
    bimFileName = '';
    famFileName = '';
    fileInput.value = '';
    badgeBed.classList.remove('loaded');
    badgeBim.classList.remove('loaded');
    badgeFam.classList.remove('loaded');
    fileName.textContent = 'No files loaded';
    inputSummary.textContent = 'No files loaded';
    runBtn.disabled = true;
    computeNote.textContent = '';
    resultsPanel.hidden = true;
    syncCLIState();
    showToast('Loaded files cleared', 'ok');
  });
}

function handleFiles(files) {
  let loadedCount = 0;
  for (const file of files) {
    const lower = file.name.toLowerCase();
    if (lower.endsWith('.bed')) {
      const reader = new FileReader();
      reader.onload = e => {
        bedBytes = new Uint8Array(e.target.result);
        bedFileName = file.name;
        badgeBed.classList.add('loaded');
        checkReady();
      };
      reader.readAsArrayBuffer(file);
      loadedCount++;
    } else if (lower.endsWith('.bim')) {
      const reader = new FileReader();
      reader.onload = e => {
        bimText = e.target.result;
        bimFileName = file.name;
        badgeBim.classList.add('loaded');
        checkReady();
      };
      reader.readAsText(file);
      loadedCount++;
    } else if (lower.endsWith('.fam')) {
      const reader = new FileReader();
      reader.onload = e => {
        famText = e.target.result;
        famFileName = file.name;
        badgeFam.classList.add('loaded');
        checkReady();
      };
      reader.readAsText(file);
      loadedCount++;
    }
  }
}

function checkReady() {
  const parts = [];
  if (bedBytes) parts.push(`BED (${(bedBytes.length / 1024).toFixed(0)} KB)`);
  if (bimText) parts.push(`BIM (${bimText.trim().split('\n').length} SNPs)`);
  if (famText) parts.push(`FAM (${famText.trim().split('\n').length} samples)`);

  if (parts.length > 0) {
    const summary = parts.join(', ');
    fileName.textContent = summary;
    inputSummary.textContent = summary;
  }

  if (bedBytes && bimText && famText) {
    runBtn.disabled = false;
    computeNote.textContent = 'Ready to execute QC';
    showToast('Genotype fileset loaded successfully', 'ok');
  } else {
    runBtn.disabled = true;
    const missing = [];
    if (!bedBytes) missing.push('.bed');
    if (!bimText) missing.push('.bim');
    if (!famText) missing.push('.fam');
    computeNote.textContent = `Awaiting ${missing.join(', ')}`;
  }
  syncCLIState();
}

for (const control of [orderSelect, roundsInput, makeBedInput]) {
  control?.addEventListener('input', syncCLIState);
  control?.addEventListener('change', syncCLIState);
}
syncCLIState();

function bytesToBase64(bytes) {
  let binary = '';
  const len = bytes.byteLength;
  for (let i = 0; i < len; i++) {
    binary += String.fromCharCode(bytes[i]);
  }
  return btoa(binary);
}

// 5. Toast Notifications
function showToast(message, type = 'ok') {
  if (!toastHost) return;
  const toast = document.createElement('div');
  toast.className = `toast ${type}`;
  toast.textContent = message;
  toastHost.appendChild(toast);
  requestAnimationFrame(() => toast.classList.add('is-visible'));
  setTimeout(() => {
    toast.classList.add('is-leaving');
    toast.classList.remove('is-visible');
    setTimeout(() => toast.remove(), 250);
  }, 2800);
}

// 6. Run QC Analysis
runBtn.addEventListener('click', async () => {
  if (!bedBytes || !bimText || !famText) return;

  runBtn.disabled = true;
  progressModal.hidden = false;
  progressFill.style.width = '15%';
  progressPercent.textContent = '15%';
  progressStage.textContent = 'Encoding inputs and preparing Wasm…';

  const mindVal = optionalThresholdValue(mindEnabledInput, mindInput, 0.1);
  const genoVal = optionalThresholdValue(genoEnabledInput, genoInput, 0.05);
  const hweVal = optionalThresholdValue(hweEnabledInput, hweInput, 0.000001);
  const mafVal = optionalThresholdValue(mafEnabledInput, mafInput, 0.05);
  const maxMafVal = optionalThresholdValue(maxMafEnabledInput, maxMafInput, 0.4);

  const validateOptional = (enabled, value, label) => {
    if (!enabled?.checked) return true;
    if (!Number.isFinite(value) || value < 0 || value > 1) {
      runBtn.disabled = false;
      progressModal.hidden = true;
      showToast(`${label} threshold must be between 0 and 1.`, 'error');
      return false;
    }
    return true;
  };
  if (!validateOptional(mindEnabledInput, mindVal, 'Mind')) return;
  if (!validateOptional(genoEnabledInput, genoVal, 'Geno')) return;
  if (!validateOptional(hweEnabledInput, hweVal, 'HWE')) return;
  if (!validateOptional(mafEnabledInput, mafVal, 'MAF')) return;
  if (!validateOptional(maxMafEnabledInput, maxMafVal, 'Max-MAF')) return;

  const payload = {
    arguments: {
      fam_text: famText,
      bim_text: bimText,
      bed_base64: bytesToBase64(bedBytes),
      mind: mindVal,
      geno: genoVal,
      hwe: hweVal,
      maf: mafVal,
      max_maf: maxMafVal,
      order: orderSelect.value,
      rounds: parseInt(roundsInput.value, 10) || 1,
      make_bed: makeBedInput.checked,
    },
    request: { run_id: 'qc-' + Date.now() },
  };

  progressFill.style.width = '45%';
  progressPercent.textContent = '45%';
  progressStage.textContent = 'Scanning genotype matrix in WebAssembly…';

  window.onGenoKitAnalysisComplete = data => {
    if (progressWatchdog) {
      clearTimeout(progressWatchdog);
      progressWatchdog = null;
    }
    progressFill.style.width = '100%';
    progressPercent.textContent = '100%';
    progressStage.textContent = 'Completed!';
    setTimeout(() => {
      progressModal.hidden = true;
      runBtn.disabled = false;
      if (data.status === 'success') {
        renderResults(data);
        persistQcArtifacts(data).catch(err => {
          showToast('QC finished but workspace save failed: ' + (err?.message || String(err)), 'error');
        });
        showToast('QC completed successfully', 'ok');
      } else {
        showToast('QC execution failed: ' + (data.error?.message || 'Unknown error'), 'error');
      }
    }, 350);
  };

  window.onPlinkitAnalysisComplete = window.onGenoKitAnalysisComplete;

  let progressWatchdog = setTimeout(() => {
    progressWatchdog = null;
    progressModal.hidden = true;
    runBtn.disabled = false;
    showToast('QC timed out in the browser Wasm path. Try CLI: genokit --bfile <prefix>', 'error');
  }, 120000);

  if (window.xmosomeBrowserExecute) {
    Promise.resolve(window.xmosomeBrowserExecute(payload)).catch(err => {
      if (progressWatchdog) {
        clearTimeout(progressWatchdog);
        progressWatchdog = null;
      }
      progressModal.hidden = true;
      runBtn.disabled = false;
      showToast('QC execution failed: ' + (err?.message || String(err)), 'error');
    });
  } else {
    if (progressWatchdog) {
      clearTimeout(progressWatchdog);
      progressWatchdog = null;
    }
    progressModal.hidden = true;
    runBtn.disabled = false;
    showToast('Browser Wasm execution bridge not found', 'error');
  }
});

// 7. Workspace save (PedigreeKit-compatible xmosomeSave / xmosomeSaveBatch bridge)
const saveJobQueue = [];
let saveJobActive = null;

function finishSaveJob(ok, message) {
  const job = saveJobActive;
  saveJobActive = null;
  if (job) {
    if (ok) job.resolve(message);
    else job.reject(new Error(message || 'Save failed'));
  } else if (message) {
    notice(message, !ok);
  }
  pumpSaveJobs();
}

function pumpSaveJobs() {
  if (saveJobActive || !saveJobQueue.length) return;
  const job = saveJobQueue.shift();
  saveJobActive = job;
  if (window.webkit?.messageHandlers?.xmosomeSave) {
    const payload = {
      kind: job.kind,
      filename: job.filename,
      content: job.content,
    };
    if (job.encoding) payload.encoding = job.encoding;
    if (job.announce) payload.announce = '1';
    if (job.open === false) payload.open = '0';
    window.webkit.messageHandlers.xmosomeSave.postMessage(payload);
    if (!job.quiet && !job.announce) notice(`Saving ${job.filename}…`);
    return;
  }
  // Browser fallback: trigger a download.
  let blob;
  if (job.encoding === 'base64') {
    const raw = atob(job.content || '');
    const arr = new Uint8Array(raw.length);
    for (let i = 0; i < raw.length; i++) arr[i] = raw.charCodeAt(i);
    blob = new Blob([arr], { type: 'application/octet-stream' });
  } else {
    const type = job.kind === 'report' ? 'text/markdown' : 'text/plain';
    blob = new Blob([job.content], { type });
  }
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = job.filename;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  if (!job.quiet) notice(`Downloaded ${job.filename}.`);
  finishSaveJob(true, `Downloaded ${job.filename}.`);
}

function saveOutput(kind, filename, content, { quiet = false, announce = false, open, encoding } = {}) {
  return new Promise((resolve, reject) => {
    saveJobQueue.push({ kind, filename, content, quiet, announce, open, encoding, resolve, reject });
    pumpSaveJobs();
  });
}

function saveBatch(files, { announceFilename = '' } = {}) {
  return new Promise((resolve, reject) => {
    if (!files.length) {
      resolve('Nothing to save');
      return;
    }
    if (!window.webkit?.messageHandlers?.xmosomeSaveBatch) {
      (async () => {
        for (const file of files) {
          await saveOutput(file.kind, file.filename, file.content, {
            quiet: true,
            open: false,
            announce: false,
            encoding: file.encoding,
          });
        }
        resolve(`Downloaded ${files.length} files.`);
      })().catch(reject);
      return;
    }
    saveJobActive = {
      resolve,
      reject,
      quiet: true,
      announce: false,
    };
    window.webkit.messageHandlers.xmosomeSaveBatch.postMessage({
      files,
      announce: announceFilename || '',
    });
  });
}

window.xmosomeSaveResult = (ok, message) => {
  if (saveJobActive) finishSaveJob(!!ok, message);
  else notice(message, !ok);
};

function buildQcSaveFiles(data) {
  const m = data.metrics || {};
  const retainedSamples = m.retained_samples ?? '—';
  const retainedVariants = m.retained_variants ?? '—';
  const noteName = 'report.md';
  const slimReport = [
    '---',
    'title: GenoKit QC',
    'tags:',
    '  - plugin',
    '  - GenoKit',
    '  - qc',
    '---',
    '',
    '# GenoKit QC',
    '',
    data.message || 'QC completed.',
    '',
    `- Retained samples: **${retainedSamples}**`,
    `- Retained variants: **${retainedVariants}**`,
    '',
    '## Quality control',
    '',
    '![QC summary](data/qc-summary.json)',
    '',
    '## Filtered fileset',
    '',
    '![Filtered FAM](data/filtered.fam)',
    '',
    '![Filtered BIM](data/filtered.bim)',
    '',
  ].join('\n');
  const files = [];
  files.push({
    kind: 'report',
    filename: noteName,
    content: slimReport,
    open: '0',
    announce: '1',
  });
  if (data.metrics) {
    files.push({
      kind: 'data',
      filename: 'qc-summary.json',
      content: JSON.stringify(data.metrics, null, 2),
      open: '0',
    });
  }
  if (data.files?.filtered_fam) {
    files.push({
      kind: 'data',
      filename: 'filtered.fam',
      content: data.files.filtered_fam,
      open: '0',
    });
  }
  if (data.files?.filtered_bim) {
    files.push({
      kind: 'data',
      filename: 'filtered.bim',
      content: data.files.filtered_bim,
      open: '0',
    });
  }
  if (data.files?.filtered_bed_base64) {
    files.push({
      kind: 'data',
      filename: 'filtered.bed',
      content: data.files.filtered_bed_base64,
      encoding: 'base64',
      open: '0',
    });
  }
  return { files, noteName };
}

async function persistQcArtifacts(data) {
  const { files, noteName } = buildQcSaveFiles(data);
  const message = await saveBatch(files, { announceFilename: noteName });
  notice(message || `Saved ${files.length} GenoKit artifacts to the workspace.`);
  return message;
}

// 8. Render Results
function renderResults(data) {
  resultsPanel.hidden = false;
  resultsPanel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });

  const m = data.metrics || {};
  statSamples.textContent = `${m.retained_samples} / ${m.initial_samples}`;
  statMind.textContent = `${m.removed_by_mind ?? 0}`;
  statVariants.textContent = `${m.retained_variants} / ${m.initial_variants}`;
  statGeno.textContent = `${m.removed_by_geno ?? 0}`;
  statHwe.textContent = `${m.removed_by_hwe ?? 0}`;
  statMaf.textContent = `${m.removed_by_maf ?? 0}`;

  // Markdown Report Preview
  if (data.report?.content) {
    reportContainer.innerHTML = formatMarkdown(data.report.content);
  } else {
    reportContainer.textContent = 'No report produced.';
  }

  // Optional local download copies (workspace save is primary, PedigreeKit-style).
  downloadButtons.innerHTML = '';
  function addDownload(filename, content, type) {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const btn = document.createElement('a');
    btn.href = url;
    btn.download = filename;
    btn.className = 'downloadBtn';
    btn.innerHTML = `<svg width="12" height="12" viewBox="0 0 16 16" fill="currentColor"><path d="M8 12l-4-4h2.5V2h3v6H12L8 12zm-6 2h12v1.5H2V14z"/></svg> ${filename}`;
    downloadButtons.appendChild(btn);
  }

  if (data.report?.content) {
    addDownload('report.md', data.report.content, 'text/markdown');
  }
  if (data.metrics) {
    addDownload('qc-summary.json', JSON.stringify(data.metrics, null, 2), 'application/json');
  }
  if (data.files?.filtered_fam) {
    addDownload('filtered.fam', data.files.filtered_fam, 'text/plain');
  }
  if (data.files?.filtered_bim) {
    addDownload('filtered.bim', data.files.filtered_bim, 'text/plain');
  }
  if (data.files?.filtered_bed_base64) {
    const raw = atob(data.files.filtered_bed_base64);
    const arr = new Uint8Array(raw.length);
    for (let i = 0; i < raw.length; i++) {
      arr[i] = raw.charCodeAt(i);
    }
    addDownload('filtered.bed', arr, 'application/octet-stream');
  }

  const saveAgain = document.createElement('button');
  saveAgain.type = 'button';
  saveAgain.className = 'downloadBtn';
  saveAgain.textContent = 'Save to workspace again';
  saveAgain.addEventListener('click', () => {
    persistQcArtifacts(data).catch(err => {
      showToast('Workspace save failed: ' + (err?.message || String(err)), 'error');
    });
  });
  downloadButtons.appendChild(saveAgain);
}

// Simple markdown formatter for clean table and heading rendering
function formatMarkdown(md) {
  let html = md
    .replace(/^### (.*$)/gim, '<h3>$1</h3>')
    .replace(/^## (.*$)/gim, '<h2>$1</h2>')
    .replace(/^# (.*$)/gim, '<h1>$1</h1>')
    .replace(/\*\*(.*?)\*\*/gim, '<strong>$1</strong>')
    .replace(/\*(.*?)\*/gim, '<em>$1</em>')
    .replace(/^---$/gim, '<hr>');

  // Markdown tables
  const lines = html.split('\n');
  let inTable = false;
  let tableHtml = '';
  const outLines = [];

  for (let line of lines) {
    line = line.trim();
    if (line.startsWith('|') && line.endsWith('|')) {
      if (line.includes('---')) continue; // separator
      const cells = line.split('|').slice(1, -1).map(c => c.trim());
      if (!inTable) {
        inTable = true;
        tableHtml = '<table><thead><tr>' + cells.map(c => `<th>${c}</th>`).join('') + '</tr></thead><tbody>';
      } else {
        tableHtml += '<tr>' + cells.map(c => `<td>${c}</td>`).join('') + '</tr>';
      }
    } else {
      if (inTable) {
        tableHtml += '</tbody></table>';
        outLines.push(tableHtml);
        inTable = false;
        tableHtml = '';
      }
      if (line.length > 0) {
        outLines.push(line.startsWith('<') ? line : `<p>${line}</p>`);
      }
    }
  }
  if (inTable) {
    tableHtml += '</tbody></table>';
    outLines.push(tableHtml);
  }
  return outLines.join('\n');
}
