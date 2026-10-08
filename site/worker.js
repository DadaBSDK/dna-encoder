// Runs error-tolerant decoding (and channel simulation) off the main thread.
import { decodeResilient, damage, ResilientError } from './resilient.js';

self.onmessage = ({ data }) => {
  try {
    if (data.op === 'damage') {
      self.postMessage({ ok: true, ...damage(data.text, data.channel) });
      return;
    }
    const out = decodeResilient(data.text);
    const bytes = out.bytes.slice();
    self.postMessage({ ok: true, bytes, verified: out.verified, meta: out.meta, report: out.report }, [bytes.buffer]);
  } catch (e) {
    self.postMessage({ ok: false, error: e instanceof ResilientError ? e.message : `decoder error: ${e.message}` });
  }
};
