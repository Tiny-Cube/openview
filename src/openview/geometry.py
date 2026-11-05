"""Geometric edits with atom identity and connectivity preserved."""
import numpy as np


def geometry_value(positions, indices):
    p = np.asarray(positions, dtype='f8')[indices]
    if len(p) == 2:
        return float(np.linalg.norm(p[1] - p[0]))
    if len(p) == 3:
        a, b = p[0] - p[1], p[2] - p[1]
        denominator = np.linalg.norm(a) * np.linalg.norm(b)
        if denominator < 1e-12:
            raise ValueError('Angle is undefined for coincident atoms.')
        return float(np.degrees(np.arccos(np.clip(a @ b / denominator, -1, 1))))
    if len(p) == 4:
        axis = p[2] - p[1]
        norm = np.linalg.norm(axis)
        if norm < 1e-12:
            raise ValueError('Dihedral axis has zero length.')
        axis /= norm
        a, b = p[0] - p[1], p[3] - p[2]
        a -= (a @ axis) * axis; b -= (b @ axis) * axis
        if min(np.linalg.norm(a), np.linalg.norm(b)) < 1e-12:
            raise ValueError('Dihedral is undefined for collinear atoms.')
        return float(np.degrees(np.arctan2(np.cross(axis, a) @ b, a @ b)))
    raise ValueError('Select 2, 3 or 4 atoms in order.')


def edited_geometry(positions, bonds, indices, target, move_group=True):
    """Fix the first side, translate/rotate the last atom or its attached group.

    A group is separated by cutting the edited bond/axis. Ring-connected groups
    cannot be separated this way and must use the explicit atom-only option.
    """
    if len(indices) not in (2, 3, 4) or len(set(indices)) != len(indices):
        raise ValueError('Select 2, 3 or 4 distinct atoms in order.')
    xyz = np.array(positions, dtype='f8', copy=True)
    current = geometry_value(xyz, indices)
    if not np.isfinite(target):
        raise ValueError('Enter a finite geometry value.')
    count = len(indices)
    if count == 2 and target <= 0:
        raise ValueError('Distance must be positive.')
    if count == 3 and not 0 < target < 180:
        raise ValueError('Angle must lie between 0 and 180 degrees.')
    if abs(target - current) < 1e-10:
        return xyz.astype(positions.dtype)
    cut = tuple(indices[:2] if count == 2 else indices[1:3])
    moving = {indices[-1]}
    if move_group:
        neighbors = [[] for _ in xyz]
        for a, b, _ in bonds:
            a, b = int(a), int(b)
            if {a, b} != set(cut):
                neighbors[a].append(b); neighbors[b].append(a)
        pending = list(moving)
        while pending:
            for other in neighbors[pending.pop()]:
                if other not in moving:
                    moving.add(other); pending.append(other)
        fixed = set(indices[:-1] if count < 4 else indices[:2])
        if moving & fixed:
            raise ValueError('This group is connected through a ring or alternate path. Disable "Move attached group" to move only the last atom.')
    selected = sorted(moving)
    if count == 2:
        direction = xyz[indices[1]] - xyz[indices[0]]
        if current < 1e-12:
            raise ValueError('Cannot choose a direction for coincident atoms.')
        xyz[selected] += direction * ((target - current) / current)
    else:
        if count == 3:
            origin = xyz[indices[1]].copy()
            axis = np.cross(xyz[indices[0]] - origin, xyz[indices[2]] - origin)
            if np.linalg.norm(axis) < 1e-12:
                # A linear angle has no unique bending plane. Choose a stable perpendicular.
                direction = xyz[indices[2]] - origin
                basis = np.eye(3)[np.argmin(np.abs(direction))]
                axis = np.cross(direction, basis)
        else:
            origin = xyz[indices[2]].copy()
            axis = xyz[indices[2]] - xyz[indices[1]]
        axis /= np.linalg.norm(axis)
        angle = np.radians(target - current)
        relative = xyz[selected] - origin
        xyz[selected] = origin + (relative * np.cos(angle) + np.cross(axis, relative) * np.sin(angle)
            + np.outer(relative @ axis, axis) * (1 - np.cos(angle)))
    return xyz.astype(positions.dtype)
