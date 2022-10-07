import wave
import contextlib
files = os.listdir('.')

for fname in files:
    with contextlib.closing(wave.open(fname,'r')) as f:
        frames = f.getnframes()
        rate = f.getframerate()
        duration = frames / float(rate)
        print(fname)
        print(duration)
