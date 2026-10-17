"""Independent, exception-based geometry and native artifact validation."""
import math
import re

from .models import matrix

NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?|[-+]?(?:nan|inf(?:inity)?)'


def finite_number(value):
    if isinstance(value, bool):
        raise ValueError('Expected finite scalar, not boolean')
    number = float(str(value).replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(number):
        raise ValueError('Expected finite scalar')
    return number


def close_matrix(first, second, *, tolerance=5.1e-6):
    return len(first) == len(second) and all(
        len(a) == len(b) == 3 and all(abs(x-y) <= tolerance for x,y in zip(a,b))
        for a,b in zip(first,second))


def read_xyz_frames(path, symbols):
    """Strict complete ordered frames. No optional chemistry packages are imported."""
    lines = path.read_text().splitlines()
    cursor, frames = 0, []
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        count = int(lines[cursor])
        if count != len(symbols) or cursor + count + 2 > len(lines):
            raise ValueError('Wrong atom count or truncated XYZ')
        positions, actual = [], []
        for line in lines[cursor+2:cursor+count+2]:
            parts = line.split()
            if len(parts) != 4:
                raise ValueError('Malformed XYZ coordinate row')
            actual.append(parts[0].capitalize())
            positions.append([finite_number(x) for x in parts[1:]])
        if tuple(actual) != symbols:
            raise ValueError('Ordered atom symbols do not match the request')
        frames.append((matrix(positions, count), lines[cursor+1]))
        cursor += count + 2
    if not frames:
        raise ValueError('Missing XYZ frames')
    return tuple(frames)


def xyz_energy(comment):
    # xTB energy label or geomeTRIC's Iteration N Energy E; never use a
    # trailing version number/gradient as the energy.
    matches = re.findall(r'(?:\benergy\s*:?\s*)(' + NUMBER + ')', comment, re.I)
    if not matches:
        raise ValueError('No frame-associated energy in XYZ comment')
    return finite_number(matches[-1])


def orca_coordinates(text, count):
    header = 'CARTESIAN COORDINATES (ANGSTROEM)'
    if header not in text:
        raise ValueError('Missing final ORCA Cartesian geometry')
    tail = text.rsplit(header, 1)[1].splitlines()
    actual, positions = [], []
    for line in tail:
        parts = line.split()
        if len(parts) == 4 and parts[0].isalpha():
            actual.append(parts[0].capitalize())
            positions.append([finite_number(x) for x in parts[1:]])
            if len(actual) == count:
                return tuple(actual), matrix(positions, count)
        elif actual:
            break
    raise ValueError('Incomplete final ORCA Cartesian geometry')


def verify_start(request, root, path):
    if request.optimizer == 'native' and request.backend == 'orca':
        text = path('trial_' + request.name + '.inp').read_text()
        match = re.search(r'\*\s*xyz\s+(-?\d+)\s+(\d+)\s*\n(.*?)\n\s*\*', text, re.I | re.S)
        if match is None:
            raise ValueError('Missing ORCA input geometry')
        settings = request.settings
        if (int(match[1]), int(match[2])) != (settings['charge'], settings['multiplicity']):
            raise ValueError('ORCA input electronic state mismatch')
        rows = [line.split() for line in match[3].splitlines() if line.strip()]
        if tuple(row[0].capitalize() for row in rows) != request.symbols:
            raise ValueError('ORCA input ordered symbols mismatch')
        coordinates = matrix([[finite_number(x) for x in row[1:]] for row in rows], len(request.symbols))
    else:
        coordinates = read_xyz_frames(path('trial_' + request.name + '.xyz'), request.symbols)[0][0]
    if not close_matrix(coordinates, request.coordinates_angstrom):
        raise ValueError('Input geometry does not match original request')
