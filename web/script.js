const API = "http://localhost:8000";

// unwrap scalar-like arrays from the R JSON (e.g. [true] -> true)
function unwrap(v){
  if (Array.isArray(v) && v.length === 1) return v[0];
  return v;
}

async function addGene(){
  const gene = document.getElementById('gene').value.trim();
  const res = await fetch(`${API}/genes`, {
    method: 'POST',
    headers: {'Content-Type':'application/json'},
    body: JSON.stringify({gene})
  });
  const data = await res.json();

  // user-friendly summary (fall back to full JSON if unexpected)
  const ok = unwrap(data.ok);
  const found = unwrap(data.found);
  const added = unwrap(data.added);
  let msg = '';

  if (ok === true) {
    if (added) msg = `Added: ${added}`;
    else if (found === false) msg = `Not found: ${gene}`;
    else msg = 'OK';
  } else {
    // show full JSON if something unexpected happened
    msg = JSON.stringify(data, null, 2);
  }

  document.getElementById('out').textContent = msg;
  refreshPlot();

  // clear the input and return focus
  const input = document.getElementById('gene');
  input.value = '';
  input.focus();

  // also keep full response in console for debugging
  console.debug('addGene response:', data);
}


async function resetGenes(){
  if (!confirm("Clear all interesting genes?")) return;
  const res = await fetch(`${API}/genes/reset`, { method: 'POST' });
  const data = await res.json();

  const cleared = unwrap(data.cleared);
  let msg = cleared ? 'All interesting genes cleared.' : JSON.stringify(data, null, 2);

  document.getElementById('out').textContent = msg;
  refreshPlot();
  console.debug('resetGenes response:', data);
}

async function uploadCSV(){
  const fileInput = document.getElementById('csvFile');
  const file = fileInput.files[0];
  
  if (!file) {
    document.getElementById('out').textContent = 'Please select a CSV file first.';
    return;
  }
  
  // show loading message
  document.getElementById('out').textContent = 'Uploading and processing CSV...';
  
  try {
    // Read the file as text
    const fileContent = await file.text();
    
    const res = await fetch(`${API}/upload`, {
      method: 'POST',
      headers: {
        'Content-Type': 'text/csv'
      },
      body: fileContent
    });
    
    const data = await res.json();
    const ok = unwrap(data.ok);
    
    if (ok === true) {
      const rows = unwrap(data.rows);
      document.getElementById('out').textContent = `CSV uploaded! Loaded ${rows} rows.`;
      refreshPlot();
    } else {
      const message = unwrap(data.message) || 'Upload failed';
      document.getElementById('out').textContent = `Error: ${message}`;
    }
    
    console.debug('uploadCSV response:', data);
    
  } catch (error) {
    document.getElementById('out').textContent = `Upload error: ${error.message}`;
    console.error('Upload error:', error);
  }
}

// Refresh plot using slider values (lfc and minlogp) as query params
function refreshPlot(){
  const lfc = encodeURIComponent(document.getElementById('lfc').value);
  const minlogp = encodeURIComponent(document.getElementById('minlogp').value);
  // cache-bust and include thresholds and whether threshold labels should be shown
  const show = (typeof window.showThresholdLabels === 'undefined') ? true : !!window.showThresholdLabels;
  document.getElementById('plot').src = `${API}/plot?lfc=${lfc}&minlogp=${minlogp}&show_threshold_labels=${show}&ts=${Date.now()}`;
}

// Update the display of converted threshold values
function updateThresholdDisplay(){
  const log2fc = parseFloat(document.getElementById('lfc').value);
  const minlogp = parseFloat(document.getElementById('minlogp').value);
  
  // Convert log2 fold change to fold change: 2^log2fc
  const foldChange = Math.pow(2, log2fc);
  
  // Convert -log10(padj) to padj: 10^(-minlogp)
  const padj = Math.pow(10, -minlogp);
  
  // Update display with appropriate formatting
  document.getElementById('fc_display').textContent = foldChange.toFixed(2);
  document.getElementById('padj_display').textContent = padj.toExponential(2);
}

// Wire sliders: update displayed numeric value and refresh the plot on change
function initSliders(){
  const lfcSlider = document.getElementById('lfc');
  const lfcNumber = document.getElementById('lfc_val');
  const minlogpSlider = document.getElementById('minlogp');
  const minlogpNumber = document.getElementById('minlogp_val');

  if (lfcSlider && lfcNumber) {
    // initialize
    lfcNumber.value = lfcSlider.value;
    // slider -> number
    lfcSlider.addEventListener('input', function(){
      lfcNumber.value = this.value;
      updateThresholdDisplay();
      refreshPlot();
    });
    // number -> slider
    lfcNumber.addEventListener('input', function(){
      lfcSlider.value = this.value;
      updateThresholdDisplay();
      refreshPlot();
    });
  }

  if (minlogpSlider && minlogpNumber) {
    // initialize
    minlogpNumber.value = minlogpSlider.value;
    // slider -> number
    minlogpSlider.addEventListener('input', function(){
      minlogpNumber.value = this.value;
      updateThresholdDisplay();
      refreshPlot();
    });
    // number -> slider (with bounds clamping)
    minlogpNumber.addEventListener('input', function(){
      const v = Math.min(Math.max(parseFloat(this.value || 0), parseFloat(minlogpSlider.min)), parseFloat(minlogpSlider.max));
      this.value = v;
      minlogpSlider.value = v;
      updateThresholdDisplay();
      refreshPlot();
    });
  }
  
  // Initialize the threshold display
  updateThresholdDisplay();
}

// initial
// default toggle state (false = hide threshold-derived labels by default)
window.showThresholdLabels = false;

function toggleThresholdLabels(){
  window.showThresholdLabels = !window.showThresholdLabels;
  const btn = document.getElementById('toggleLabelsBtn');
  if (btn) btn.textContent = `Threshold labels: ${window.showThresholdLabels ? 'On' : 'Off'}`;
  refreshPlot();
}

// wire the toggle button if present
document.addEventListener('DOMContentLoaded', function(){
  const btn = document.getElementById('toggleLabelsBtn');
  if (btn) {
    btn.addEventListener('click', toggleThresholdLabels);
    btn.textContent = `Threshold labels: ${window.showThresholdLabels ? 'On' : 'Off'}`;
  }
});

initSliders();
refreshPlot();

// call addGene() when Enter is pressed in the input
document.getElementById('gene').addEventListener('keydown', function(e){
  if (e.key === 'Enter') {
    e.preventDefault();
    addGene();
  }
});

(function wireThresholdInputs(){
  const lfcSlider = document.getElementById('lfc');
  const lfcNumber = document.getElementById('lfc_val');
  const minlogpSlider = document.getElementById('minlogp');
  const minlogpNumber = document.getElementById('minlogp_val');

  if (lfcSlider && lfcNumber) {
    // initialize
    lfcNumber.value = lfcSlider.value;
    // slider -> number
    lfcSlider.addEventListener('input', () => {
      lfcNumber.value = lfcSlider.value;
      refreshPlot();
    });
    // number -> slider
    lfcNumber.addEventListener('input', () => {
      lfcSlider.value = lfcNumber.value;
      refreshPlot();
    });
  }

  if (minlogpSlider && minlogpNumber) {
    minlogpNumber.value = minlogpSlider.value;
    minlogpSlider.addEventListener('input', () => {
      minlogpNumber.value = minlogpSlider.value;
      refreshPlot();
    });
    minlogpNumber.addEventListener('input', () => {
      // clamp to slider bounds
      const v = Math.min(Math.max(parseFloat(minlogpNumber.value || 0), parseFloat(minlogpSlider.min)), parseFloat(minlogpSlider.max));
      minlogpNumber.value = v;
      minlogpSlider.value = v;
      refreshPlot();
    });
  }
})();