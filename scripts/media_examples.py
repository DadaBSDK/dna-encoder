"""Create valid synthetic media fixtures and verify streaming DNA recovery.

Requires ffmpeg for fixture generation only; dnastore itself does not require it.
"""
from pathlib import Path
import json
import subprocess

from dnastore.stream import encode_file, decode_file


def main():
    root = Path('results/media_stream_examples')
    root.mkdir(parents=True, exist_ok=False)
    commands = {
        'image.png': ['-f', 'lavfi', '-i', 'testsrc2=size=320x240', '-frames:v', '1'],
        'animation.gif': ['-f', 'lavfi', '-i', 'testsrc2=size=96x64:rate=6', '-t', '2'],
        'video.mp4': ['-f', 'lavfi', '-i', 'testsrc2=size=320x240:rate=24', '-t', '8', '-c:v', 'libx264', '-crf', '18', '-pix_fmt', 'yuv420p'],
    }
    results = []
    for name, args in commands.items():
        original = root / name
        subprocess.run(['ffmpeg', '-nostdin', '-v', 'error', *args, str(original)], check=True)
        archive = root / (name + '.dna')
        recovered = root / ('recovered-' + name)
        record = encode_file(original, archive)
        assert decode_file(archive, recovered) == record
        assert original.read_bytes() == recovered.read_bytes()
        # Decode media using ffmpeg as an additional format validity check.
        subprocess.run(['ffmpeg', '-nostdin', '-v', 'error', '-i', str(recovered), '-f', 'null', '-'], check=True)
        results.append({**record, 'exact_recovery': True, 'media_decode': 'passed'})
    (root / 'results.json').write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
