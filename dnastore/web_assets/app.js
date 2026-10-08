'use strict';
const $ = id => document.getElementById(id);
let selectedFile = null, token = '', current = null, jobs = new Map(), busy = false, previewURL = null, started = 0;
const num = n => Number(n || 0).toLocaleString();
const bytes = n => n < 1024 ? `${n} B` : n < 1048576 ? `${(n/1024).toFixed(1)} KB` : `${(n/1048576).toFixed(1)} MB`;
const mode = () => document.querySelector('input[name=mode]:checked').value;
const show = (id, yes) => {$(id).hidden = !yes;};
function error(message) {$('form-error').textContent = message; show('form-error', !!message);}
function setFile(file) {
  selectedFile = file; error(''); $('run').disabled = busy || !file || !token;
  $('file-title').textContent = file ? file.name : 'Drop a file to begin';
  $('file-subtitle').textContent = file ? `${bytes(file.size)} · ready to upload` : 'Images, videos, GIFs, documents — any file';
  $('file-meta').textContent = file ? `${bytes(file.size)} input${mode() === 'stream' && $('action').value !== 'decode' ? ` → ${num(file.size*4)} DNA nucleotides` : ''}` : '';
  show('file-meta', !!file); $('media-preview').replaceChildren(); show('media-preview', false);
  if (previewURL) URL.revokeObjectURL(previewURL);
  previewURL = null;
  if (file && /^(image|video|audio)\//.test(file.type) && !file.type.includes('svg')) {
    const element = document.createElement(file.type.startsWith('image') ? 'img' : file.type.startsWith('video') ? 'video' : 'audio');
    previewURL = URL.createObjectURL(file); element.src = previewURL;
    if (element.tagName !== 'IMG') {element.controls = true; element.preload = 'metadata';} else element.alt = 'Selected file preview';
    $('media-preview').append(element); show('media-preview', true);
  }
}
function settings() {
  const oligo = mode() === 'oligo';
  show('oligo-settings', oligo);
  $('mode-note').textContent = oligo ? 'Indexed oligos, error correction, and controlled sequencing experiments. Best for small test files.' : 'Exact digital storage with SHA-256 verification. Uses about 4× the original disk space.';
  const roundtrip = $('action').value === 'roundtrip';
  $('simulate').disabled = !roundtrip;
  show('channel-settings', oligo && roundtrip && $('simulate').checked);
  $('codec').disabled = $('action').value === 'decode';
  $('parity').disabled = $('action').value === 'decode' || $('codec').value === 'fountain';
  $('file').accept = $('action').value === 'decode' ? (oligo ? '.fasta,.fa,.fastq,.fq,.gz' : '.dna') : '';
  if (selectedFile) $('file-meta').textContent = `${bytes(selectedFile.size)} input${!oligo && $('action').value !== 'decode' ? ` → ${num(selectedFile.size*4)} DNA nucleotides` : ''}`;
}
$('file').addEventListener('change', e => setFile(e.target.files[0] || null));
for (const event of ['dragenter','dragover']) $('dropzone').addEventListener(event, e => {e.preventDefault(); $('dropzone').classList.add('drag');});
for (const event of ['dragleave','drop']) $('dropzone').addEventListener(event, e => {e.preventDefault(); $('dropzone').classList.remove('drag'); if(event === 'drop') setFile(e.dataTransfer.files[0] || null);});
$('sample').onclick = () => {$('action').value = 'roundtrip'; settings(); setFile(new File(['DNA storage turns bytes into bases.\nImages, videos, GIFs and documents share the same journey.\nRecover every byte, then verify it.\n'], 'hello-dna.txt', {type:'text/plain'}));};
for (const input of document.querySelectorAll('input[name=mode], #action, #simulate, #codec')) input.addEventListener('change', settings);
for (const b of document.querySelectorAll('[data-preset]')) b.onclick = () => {const p = b.dataset.preset; for(const id of ['sub','ins','del']) $(id).value = p === 'clean' ? 0 : p === 'moderate' ? 1 : 3; $('depth').value = p === 'stress' ? 2 : 12; $('dropout').value = 0;};
function options() {
  return {mode:mode(), action:$('action').value, codec:$('codec').value, oligo_len:Number($('oligo-len').value), parity:Number($('parity').value), seed:Number($('seed').value), strength:$('strength').value, simulate:$('simulate').checked && $('action').value === 'roundtrip', channel:{p_sub:Number($('sub').value)/100,p_ins:Number($('ins').value)/100,p_del:Number($('del').value)/100,dropout:Number($('dropout').value)/100,mean_coverage:Number($('depth').value)}};
}
function upload(file, opts) {
  return new Promise((resolve,reject) => {
    const xhr = new XMLHttpRequest(); xhr.open('POST','/api/run');
    xhr.setRequestHeader('X-Workbench-Token',token); xhr.setRequestHeader('X-Filename',encodeURIComponent(file.name)); xhr.setRequestHeader('X-Options',encodeURIComponent(JSON.stringify(opts)));
    xhr.upload.onprogress = e => {if(e.lengthComputable) $('upload-progress').value = e.loaded/e.total*100;};
    xhr.upload.onload = () => {$('phase').textContent='Upload complete · preparing experiment'; $('upload-progress').removeAttribute('value');};
    xhr.onload = () => {try {const body=JSON.parse(xhr.responseText); xhr.status===202 ? resolve(body) : reject(new Error(body.error || 'Upload failed'));} catch {reject(new Error('Server returned an invalid response'));}};
    xhr.onerror=()=>reject(new Error('Connection lost. Check that the local server is running.')); xhr.send(file);
  });
}
function history() {
  $('run-count').textContent=jobs.size; if(!jobs.size)return;
  $('history-list').replaceChildren();
  for(const job of [...jobs.values()].reverse()) {
    const row=document.createElement('button'); row.className='history-row'; row.type='button';
    for(const [tag,text] of [['strong',job.name],['span',job.options.mode==='stream'?'Streaming':'Oligo lab'],['span',new Date(job.created*1000).toLocaleTimeString([], {hour:'2-digit',minute:'2-digit'})],['span',job.status==='done'?(job.result.verified===false?'Not recovered':'Complete'):job.status]]) {const n=document.createElement(tag); n.textContent=text; row.append(n);}
    row.onclick=()=>{if(busy){error('Please wait for the active experiment to finish before switching runs.');return;} current=job.id; render(job); $('experiment').scrollIntoView();}; $('history-list').append(row);
  }
}
function render(job) {
  jobs.set(job.id,job); history();
  if(job.id!==current)return;
  show('empty',false); show('result-error',job.status==='error');
  const active=['running','queued'].includes(job.status); show('progress',active); show('result-content',job.status==='done');
  $('status-badge').textContent=active?'RUNNING':job.status==='error'?'ERROR':'COMPLETE'; $('phase').textContent=job.phase;
  if(job.status==='error'){$('result-error').textContent=job.error;return;} if(job.status!=='done')return;
  const r=job.result, m=r.metrics||{}, failed=r.verified===false;
  $('outcome').className='outcome'+(failed?' failed':'');
  $('outcome').textContent=failed?'Recovery failed · no verified output file':r.exact_match?'✓ Exact recovery · original and recovered SHA-256 match':r.verified?'✓ File integrity verified · recovery complete':'✓ Encoding complete · archive ready to download';
  $('metrics').replaceChildren();
  for(const [label,value,note] of [['INPUT SIZE',bytes(r.input_bytes),r.input_name],['DNA NUCLEOTIDES',m.nucleotides!==undefined?num(m.nucleotides):m.total_nt!==undefined?num(m.total_nt):'—',m.n_oligos?`${num(m.n_oligos)} oligos`:'Digital archive'],['PROCESSING TIME',`${r.elapsed_seconds}s`,'Excludes upload time']]) {
    const d=document.createElement('div');d.className='metric'; for(const [tag,text] of [['label',label],['strong',value],['small',note]]){const n=document.createElement(tag);n.textContent=text;if(tag==='small'){n.title=text;}d.append(n);} $('metrics').append(d);
  }
  $('sequence').replaceChildren(); for(const base of r.sequence||''){const n=document.createElement('span');n.className=base.toLowerCase();n.textContent=base;$('sequence').append(n);}
  const seq=r.sequence||'', gc=seq?100*[...seq].filter(c=>'GC'.includes(c)).length/seq.length:0;
  $('sequence-note').textContent=seq?`Preview GC content: ${gc.toFixed(1)}% · ${num(seq.length)} bases shown. This is a sample, not a whole-pool statistic.`:'No DNA sequence available.';
  $('checks').replaceChildren();
  const notes=[];
  if(r.exact_match!==undefined) notes.push([r.exact_match?'Original and recovered files are byte-identical':'Original and recovered checksums differ',r.recovered_sha256]);
  else if(r.recovered_sha256) notes.push(['Recovered SHA-256',r.recovered_sha256]);
  if(r.screening)notes.push([`${r.screening.violating_oligos} / ${r.screening.oligos} oligos violate screening constraints`,r.research_mode]);
  if(r.channel)notes.push([`${num(r.channel.reads)} simulated reads · ${num(r.channel.dropped_oligos)} oligos dropped`,`${num(r.channel.substitutions)} substitutions · ${num(r.channel.insertions)} insertions · ${num(r.channel.deletions)} deletions`]);
  if(r.recovery_failure)notes.push([r.recovery_failure,'']);
  for(const [text,detail] of notes){const p=document.createElement('p');p.textContent=text;const code=document.createElement('code');code.textContent=detail;p.append(code);$('checks').append(p);}
  $('downloads').replaceChildren(); for(const a of job.artifacts){const link=document.createElement('a');link.href=a.url;link.download=a.name;link.textContent=`↓ ${a.label} · ${bytes(a.bytes)}`;$('downloads').append(link);}
  $('json-report').textContent=JSON.stringify(r,null,2);
}
async function poll(id) {
  while(true){const response=await fetch(`/api/jobs/${id}`);if(!response.ok)throw new Error('Could not retrieve the experiment');const job=await response.json();render(job);if(['done','error'].includes(job.status))return;await new Promise(resolve=>setTimeout(resolve,800));}
}
$('experiment-form').onsubmit=async e=>{
  e.preventDefault();if(!selectedFile||busy)return;error('');busy=true;current=null;started=Date.now();$('run').disabled=true;$('run').textContent='Experiment running…';show('empty',false);show('result-content',false);show('result-error',false);show('progress',true);$('phase').textContent='Uploading your file';$('status-badge').textContent='UPLOADING';$('upload-progress').value=0;
  try{const job=await upload(selectedFile,options());current=job.id;await poll(job.id);}catch(err){error(err.message);show('progress',false);$('status-badge').textContent='DISCONNECTED';}finally{busy=false;$('run').textContent='Run experiment ↗';$('run').disabled=!selectedFile||!token;}
};
setInterval(()=>{if(busy)$('elapsed').textContent=`${Math.floor((Date.now()-started)/1000)}s`;},1000);
async function init(){try{const response=await fetch('/api/session');if(!response.ok)throw new Error('Local server unavailable');const data=await response.json();token=data.token;for(const j of data.jobs)jobs.set(j.id,j);history();$('connection').textContent='● Local server connected';$('run').disabled=!selectedFile;const running=data.jobs.find(j=>['running','queued'].includes(j.status));if(running){current=running.id;busy=true;started=Date.now();$('run').disabled=true;try{await poll(running.id);}finally{busy=false;$('run').disabled=!selectedFile;}}}catch(e){$('connection').textContent='Server disconnected';error(e.message);}}
settings();init();
