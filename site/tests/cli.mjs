// Node CLI used by tests/test_site_stream.py: encode/decode with the browser implementation.
import { readFileSync, writeFileSync } from 'node:fs';
import { encodeBlob, decodeBlob, StreamError } from '../stream.js';

const [cmd, src, dst, name] = process.argv.slice(2);
const blob = new Blob([readFileSync(src)]);
try {
  if (cmd === 'encode') {
    const { archive } = await encodeBlob(blob, name);
    writeFileSync(dst, Buffer.from(await archive.arrayBuffer()));
  } else {
    const { file, result } = await decodeBlob(blob);
    writeFileSync(dst, Buffer.from(await file.arrayBuffer()));
    process.stdout.write(JSON.stringify(result));
  }
} catch (e) {
  if (!(e instanceof StreamError)) throw e;
  process.stderr.write(e.message);
  process.exit(3);
}
