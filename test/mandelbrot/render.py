#!/usr/bin/env python3
"""Render the mandelbrot test's .vtu output to PNG images.

Reads the three files the test writes:
    julia-physics.vtu      escape counts on the base mesh
    julia-sensitivity.vtu  log10 |dJ/dz0| on the base mesh
    julia-refined.vtu      escape counts on the refined mesh
and produces one PNG for each, plus julia-flagged.png showing the
quarter of cells the adjoint flagged for refinement (recomputed here
with the same bisection rule the test uses).

Self-contained: numpy + PIL only, no vtk and no matplotlib.
"""

import re
import numpy as np
from PIL import Image

SIZE = 900  # output image width/height in pixels


def read_vtu(path):
    """Return (points, triangles, cell_values) from an ascii .vtu file."""
    text = open(path).read()

    def block(name, after=None):
        s = text if after is None else text[text.index(after):]
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                      s, re.S) if name else None
        return m.group(1)

    # points: the first DataArray inside <Points>
    pts_xml = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', text, re.S).group(1)
    pts = np.fromstring(pts_xml, sep=' ').reshape(-1, 3)[:, :2]

    conn_xml = re.search(r'<DataArray[^>]*Name="connectivity"[^>]*>(.*?)</DataArray>',
                         text, re.S).group(1)
    conn = np.fromstring(conn_xml, sep=' ', dtype=int)

    offs_xml = re.search(r'<DataArray[^>]*Name="offsets"[^>]*>(.*?)</DataArray>',
                         text, re.S).group(1)
    offs = np.fromstring(offs_xml, sep=' ', dtype=int)

    # cells are triangles in these files: offsets step by 3
    tris = conn.reshape(-1, 3)

    # CellData holds "volume" first, then the named field - take the last
    cd_xml = re.search(r'<CellData>(.*?)</CellData>', text, re.S).group(1)
    arrays = re.findall(r'<DataArray[^>]*>(.*?)</DataArray>', cd_xml, re.S)
    vals = np.fromstring(arrays[-1], sep=' ')

    return pts, tris, vals


def ramp(t, light, dark):
    """Single-hue sequential ramp, light -> dark, gamma-aware."""
    lo = (np.array(light, float) / 255.0) ** 2.2
    hi = (np.array(dark, float) / 255.0) ** 2.2
    c = lo[None, :] * (1.0 - t[:, None]) + hi[None, :] * t[:, None]
    return (255.0 * c ** (1 / 2.2)).astype(np.uint8)


def rasterize(pts, tris, rgb, size=SIZE):
    """Fill each triangle with its color; white background."""
    x, y = pts[:, 0], pts[:, 1]
    sx = (size - 20) / (x.max() - x.min())
    sy = (size - 20) / (y.max() - y.min())
    s = min(sx, sy)
    px = 10 + (x - x.min()) * s
    py = size - 10 - (y - y.min()) * s   # flip: image rows run downward

    img = np.full((size, size, 3), 255, np.uint8)
    for k, (a, b, c) in enumerate(tris):
        xs = np.array([px[a], px[b], px[c]])
        ys = np.array([py[a], py[b], py[c]])
        x0, x1 = int(xs.min()), int(np.ceil(xs.max()))
        y0, y1 = int(ys.min()), int(np.ceil(ys.max()))
        gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + 0.5,
                             np.arange(y0, y1 + 1) + 0.5)
        # barycentric sign tests against the three edges
        d = (xs[1] - xs[0]) * (ys[2] - ys[0]) - (xs[2] - xs[0]) * (ys[1] - ys[0])
        if d == 0:
            continue
        w0 = ((xs[1] - gx) * (ys[2] - gy) - (xs[2] - gx) * (ys[1] - gy)) / d
        w1 = ((xs[2] - gx) * (ys[0] - gy) - (xs[0] - gx) * (ys[2] - gy)) / d
        w2 = 1.0 - w0 - w1
        mask = (w0 >= -1e-9) & (w1 >= -1e-9) & (w2 >= -1e-9)
        sub = img[y0:y1 + 1, x0:x1 + 1]
        sub[mask] = rgb[k]
    return Image.fromarray(img)


def normalized(vals):
    lo, hi = vals.min(), vals.max()
    return (vals - lo) / (hi - lo) if hi > lo else np.zeros_like(vals)


def flag_top_quarter(vals):
    """The test's rule: bisect a threshold until a quarter carries it."""
    lo, hi = vals.min(), vals.max()
    n = len(vals)
    for _ in range(60):
        tau = 0.5 * (lo + hi)
        if np.count_nonzero(vals >= tau) > n // 4:
            lo = tau
        else:
            hi = tau
    return vals >= tau


def rasterize_zoom(pts, tris, rgb, win, size=700):
    """Fill triangles inside the window `win` = (x0, x1, y0, y1) and draw
    their edges, so the mesh itself is visible."""
    from PIL import ImageDraw
    x0, x1, y0, y1 = win
    cent = pts[tris].mean(axis=1)
    keep = np.nonzero((cent[:, 0] >= x0) & (cent[:, 0] <= x1) &
                      (cent[:, 1] >= y0) & (cent[:, 1] <= y1))[0]
    s = size / (x1 - x0)
    px = (pts[:, 0] - x0) * s
    py = size - (pts[:, 1] - y0) * s

    img = np.full((size, size, 3), 255, np.uint8)
    for k in keep:
        a, b, c = tris[k]
        xs = np.array([px[a], px[b], px[c]])
        ys = np.array([py[a], py[b], py[c]])
        bx0, bx1 = max(0, int(xs.min())), min(size - 1, int(np.ceil(xs.max())))
        by0, by1 = max(0, int(ys.min())), min(size - 1, int(np.ceil(ys.max())))
        if bx1 < bx0 or by1 < by0:
            continue
        gx, gy = np.meshgrid(np.arange(bx0, bx1 + 1) + 0.5,
                             np.arange(by0, by1 + 1) + 0.5)
        d = (xs[1] - xs[0]) * (ys[2] - ys[0]) - (xs[2] - xs[0]) * (ys[1] - ys[0])
        if d == 0:
            continue
        w0 = ((xs[1] - gx) * (ys[2] - gy) - (xs[2] - gx) * (ys[1] - gy)) / d
        w1 = ((xs[2] - gx) * (ys[0] - gy) - (xs[0] - gx) * (ys[2] - gy)) / d
        w2 = 1.0 - w0 - w1
        mask = (w0 >= -1e-9) & (w1 >= -1e-9) & (w2 >= -1e-9)
        sub = img[by0:by1 + 1, bx0:bx1 + 1]
        sub[mask] = rgb[k]

    im = Image.fromarray(img)
    draw = ImageDraw.Draw(im)
    for k in keep:
        a, b, c = tris[k]
        poly = [(px[a], py[a]), (px[b], py[b]), (px[c], py[c]), (px[a], py[a])]
        draw.line(poly, fill=(110, 110, 110), width=1)
    return im


def zoom_panel(win, blues):
    """Side-by-side crop of the base and refined paintings, edges drawn,
    both colored on the shared escape-count range."""
    from PIL import ImageDraw
    pts_b, tris_b, esc_b = read_vtu('julia-physics.vtu')
    pts_f, tris_f, esc_f = read_vtu('julia-refined.vtu')
    lo = min(esc_b.min(), esc_f.min())
    hi = max(esc_b.max(), esc_f.max())
    tb = (esc_b - lo) / (hi - lo)
    tf = (esc_f - lo) / (hi - lo)

    left = rasterize_zoom(pts_b, tris_b, ramp(tb, *blues), win)
    right = rasterize_zoom(pts_f, tris_f, ramp(tf, *blues), win)

    w, h, head = left.width, left.height, 26
    out = Image.new('RGB', (2 * w + 12, h + head), (255, 255, 255))
    out.paste(left, (0, head))
    out.paste(right, (w + 12, head))
    draw = ImageDraw.Draw(out)
    draw.text((6, 6), "base mesh (3700 cells)", fill=(60, 60, 60))
    draw.text((w + 18, 6), "refined, selectively marched (14800 cells)", fill=(60, 60, 60))
    out.save('julia-zoom.png')


def main():
    blues = ((222, 235, 247), (8, 48, 107))     # escape counts: light -> dark blue
    orang = ((254, 230, 206), (127, 39, 4))     # sensitivity:   light -> dark orange

    pts, tris, esc = read_vtu('julia-physics.vtu')
    rasterize(pts, tris, ramp(normalized(esc), *blues)).save('julia-physics.png')

    # one sensitivity picture per refinement cycle
    import glob
    for path in sorted(glob.glob('julia-sensitivity-*.vtu')):
        pts, tris, sens = read_vtu(path)
        out = path.replace('.vtu', '.png')
        rasterize(pts, tris, ramp(normalized(sens), *orang)).save(out)

    pts, tris, sens = read_vtu('julia-sensitivity-1.vtu')
    flags = flag_top_quarter(sens)
    rgb = np.where(flags[:, None], np.array([[178, 24, 43]]), np.array([[232, 232, 232]]))
    rasterize(pts, tris, rgb.astype(np.uint8)).save('julia-flagged.png')
    print(f"flagged {np.count_nonzero(flags)} of {len(flags)} base cells")

    pts, tris, esc = read_vtu('julia-refined.vtu')
    rasterize(pts, tris, ramp(normalized(esc), *blues)).save('julia-refined.png')

    # a crop straddling the set's boundary, centred on the flagged cell
    # nearest the middle of the canvas, edges drawn - where the selective
    # marching is visible triangle by triangle
    pts_b, tris_b, _ = read_vtu('julia-physics.vtu')
    cent = pts_b[tris_b].mean(axis=1)
    fc = cent[flags]
    mid = fc[np.argmin((fc[:, 0] - 0.55) ** 2 + (fc[:, 1] - 0.5) ** 2)]
    r = 0.14
    zoom_panel((mid[0] - r, mid[0] + r, mid[1] - r, mid[1] + r), blues)

    print("wrote julia-physics.png julia-sensitivity-<k>.png julia-flagged.png "
          "julia-refined.png julia-zoom.png")


if __name__ == '__main__':
    main()
