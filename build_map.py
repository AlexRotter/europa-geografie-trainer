"""Erzeugt mapdata.js (Europa-Karte) aus Natural-Earth-TopoJSON (world-atlas, 1:50 Mio.)."""
import json, math, sys

SRC = sys.argv[1]
OUT = sys.argv[2]

NAMES = {  # Natural-Earth-Name -> Kürzel
    'Iceland': 'ISL', 'Norway': 'NOR', 'Sweden': 'SWE', 'Finland': 'FIN', 'Denmark': 'DNK',
    'Estonia': 'EST', 'Latvia': 'LVA', 'Lithuania': 'LTU', 'United Kingdom': 'GBR', 'Ireland': 'IRL',
    'France': 'FRA', 'Belgium': 'BEL', 'Netherlands': 'NLD', 'Luxembourg': 'LUX', 'Germany': 'DEU',
    'Poland': 'POL', 'Czechia': 'CZE', 'Slovakia': 'SVK', 'Austria': 'AUT', 'Switzerland': 'CHE',
    'Hungary': 'HUN', 'Portugal': 'PRT', 'Spain': 'ESP', 'Italy': 'ITA', 'Greece': 'GRC', 'Malta': 'MLT',
    'Slovenia': 'SVN', 'Croatia': 'HRV', 'Bosnia and Herz.': 'BIH', 'Serbia': 'SRB', 'Montenegro': 'MNE',
    'Kosovo': 'XKX', 'Macedonia': 'MKD', 'North Macedonia': 'MKD', 'Albania': 'ALB', 'Romania': 'ROU',
    'Bulgaria': 'BGR', 'Russia': 'RUS', 'Belarus': 'BLR', 'Ukraine': 'UKR', 'Moldova': 'MDA', 'Turkey': 'TUR',
}

LON0, LAT0 = math.radians(15), math.radians(53)


def proj(lon, lat):
    l, p = math.radians(lon), math.radians(lat)
    k = math.sqrt(2 / (1 + math.sin(LAT0) * math.sin(p) + math.cos(LAT0) * math.cos(p) * math.cos(l - LON0)))
    x = k * math.cos(p) * math.sin(l - LON0)
    y = k * (math.cos(LAT0) * math.sin(p) - math.sin(LAT0) * math.cos(p) * math.cos(l - LON0))
    return x, -y


# Bildausschnitt aus Eckpunkten bestimmen
ext = [proj(*c) for c in [(-25, 64), (-24, 66.5), (-10, 36), (56, 41), (56, 47), (30, 71.5), (25, 34.5), (45, 66)]]
minx = min(e[0] for e in ext); maxx = max(e[0] for e in ext)
miny = min(e[1] for e in ext); maxy = max(e[1] for e in ext)
W = 1000
S = W / (maxx - minx)
H = round((maxy - miny) * S)


def px(lon, lat):
    x, y = proj(lon, lat)
    return (x - minx) * S, (y - miny) * S


topo = json.load(open(SRC))
sc, tr = topo['transform']['scale'], topo['transform']['translate']
arcs = []
for a in topo['arcs']:
    x = y = 0; pts = []
    for dx, dy in a:
        x += dx; y += dy
        pts.append((x * sc[0] + tr[0], y * sc[1] + tr[1]))
    arcs.append(pts)


def arc(i):
    return arcs[i] if i >= 0 else arcs[~i][::-1]


def ring(idx):
    pts = []
    for i in idx:
        a = arc(i)
        pts.extend(a if not pts else a[1:])
    return pts


def dp(pts, eps):
    if len(pts) < 3:
        return pts
    keep = [False] * len(pts); keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        ax, ay = pts[a]; bx, by = pts[b]
        dx, dy = bx - ax, by - ay; L = math.hypot(dx, dy) or 1e-9
        best, bi = -1, -1
        for i in range(a + 1, b):
            d = abs(dy * (pts[i][0] - ax) - dx * (pts[i][1] - ay)) / L
            if d > best:
                best, bi = d, i
        if best > eps:
            keep[bi] = True
            stack += [(a, bi), (bi, b)]
    return [p for p, k in zip(pts, keep) if k]


M = 40  # Rand, auf den außerhalb liegende Punkte geklemmt werden


def ring_path(r):
    pts = [px(lo, la) for lo, la in r]
    if all(x < -M or x > W + M or y < -M or y > H + M for x, y in pts):
        if not any(-M <= x <= W + M and -M <= y <= H + M for x, y in pts):
            # komplett außerhalb (z. B. Französisch-Guayana) – prüfen per Bounding-Box
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            if max(xs) < -M or min(xs) > W + M or max(ys) < -M or min(ys) > H + M:
                return ''
    pts = [(min(max(x, -M), W + M), min(max(y, -M), H + M)) for x, y in pts]
    # geschlossener Ring: am entferntesten Punkt teilen, sonst kollabiert Douglas-Peucker
    far = max(range(len(pts)), key=lambda i: math.hypot(pts[i][0] - pts[0][0], pts[i][1] - pts[0][1]))
    pts = dp(pts[:far + 1], 0.35)[:-1] + dp(pts[far:], 0.35)
    out = []
    for x, y in pts:
        s = f'{x:.1f},{y:.1f}'
        if not out or out[-1] != s:
            out.append(s)
    if len(out) < 3:
        return ''
    xs = [float(s.split(',')[0]) for s in out]; ys = [float(s.split(',')[1]) for s in out]
    if max(xs) - min(xs) < 0.8 and max(ys) - min(ys) < 0.8:
        return ''
    return 'M' + 'L'.join(out) + 'Z'


def geom_path(g):
    if g['type'] == 'Polygon':
        polys = [g['arcs']]
    elif g['type'] == 'MultiPolygon':
        polys = g['arcs']
    else:
        return ''
    return ''.join(ring_path(ring(r)) for poly in polys for r in poly)


countries, other = {}, []
for g in topo['objects']['countries']['geometries']:
    name = g['properties']['name']
    d = geom_path(g)
    if not d:
        continue
    code = NAMES.get(name)
    if code:
        countries[code] = d
    else:
        other.append(d)

missing = sorted(set(NAMES.values()) - set(countries))
print('fehlend:', missing, 'H =', H, file=sys.stderr)

SEAS = {  # Buchstabe wie auf dem Westermann-Arbeitsblatt
    'A': [(56.3, 19.0)], 'B': [(56.0, 3.0)], 'C': [(69.0, 2.0)], 'D': [(45.0, -12.0)],
    'E': [(38.0, 4.5), (35.4, 18.6)], 'F': [(40.0, 12.0)], 'G': [(43.0, 15.6)], 'H': [(38.0, 18.6)],
    'I': [(38.7, 25.0)], 'J': [(43.3, 34.0)], 'K': [(42.0, 50.6)],
}
seas = {k: [[round(c, 1) for c in px(lon, lat)] for lat, lon in v] for k, v in SEAS.items()}

CITIES = {  # Hauptstädte (lat, lon)
    'ISL': (64.15, -21.94), 'NOR': (59.91, 10.75), 'SWE': (59.33, 18.07), 'FIN': (60.17, 24.94),
    'DNK': (55.68, 12.57), 'EST': (59.44, 24.75), 'LVA': (56.95, 24.11), 'LTU': (54.69, 25.28),
    'GBR': (51.51, -0.13), 'IRL': (53.35, -6.26), 'FRA': (48.86, 2.35), 'BEL': (50.85, 4.35),
    'NLD': (52.37, 4.90), 'LUX': (49.61, 6.13), 'DEU': (52.52, 13.40), 'POL': (52.23, 21.01),
    'CZE': (50.08, 14.44), 'SVK': (48.15, 17.11), 'AUT': (48.21, 16.37), 'CHE': (46.95, 7.45),
    'HUN': (47.50, 19.04), 'PRT': (38.72, -9.14), 'ESP': (40.42, -3.70), 'ITA': (41.90, 12.50),
    'GRC': (37.98, 23.73), 'MLT': (35.90, 14.51), 'SVN': (46.06, 14.51), 'HRV': (45.81, 15.98),
    'BIH': (43.86, 18.41), 'SRB': (44.79, 20.45), 'MNE': (42.44, 19.26), 'XKX': (42.66, 21.17),
    'MKD': (42.00, 21.43), 'ALB': (41.33, 19.82), 'ROU': (44.43, 26.10), 'BGR': (42.70, 23.32),
    'RUS': (55.76, 37.62), 'BLR': (53.90, 27.57), 'UKR': (50.45, 30.52), 'MDA': (47.01, 28.86),
    'TUR': (39.93, 32.86),
}
cities = {k: [round(c, 1) for c in px(lon, lat)] for k, (lat, lon) in CITIES.items()}

with open(OUT, 'w') as f:
    f.write('/* Europa-Karte: Natural Earth 1:50 Mio. (gemeinfrei) über world-atlas 2.0.2,\n'
            '   Lambert-Projektion (flächentreu), erzeugt mit build_map.py */\n')
    f.write('const MAP_W = %d, MAP_H = %d;\n' % (W, H))
    f.write('const MAP_COUNTRIES = ' + json.dumps(countries, separators=(',', ':')) + ';\n')
    f.write('const MAP_OTHER = ' + json.dumps(other, separators=(',', ':')) + ';\n')
    f.write('const MAP_SEAS = ' + json.dumps(seas, separators=(',', ':')) + ';\n')
    f.write('const MAP_CITIES = ' + json.dumps(cities, separators=(',', ':')) + ';\n')
