#!/usr/bin/env python3
"""
SnapGene .dna reader/writer for ppigFinder.

This module is intentionally independent from the PyQt GUI.
"""

from pathlib import Path
from datetime import datetime
import struct as _struct
import io as _io


def _snapgene_read_packets(data: bytes):
    """Iterate TLV packets from a SnapGene .dna file."""
    buf = _io.BytesIO(data)
    cookie = buf.read(1)
    if cookie != b'\x09':
        raise ValueError(
            "File not recognised as SnapGene (.dna).\n"
            "Check whether it is corrupted or an unsupported old version."
        )
    while True:
        hdr = buf.read(5)
        if len(hdr) < 5:
            break
        pkt_type = hdr[0]
        pkt_len  = _struct.unpack('>I', hdr[1:5])[0]
        payload  = buf.read(pkt_len)
        yield pkt_type, payload


def parse_snapgene_dna(filepath: str) -> dict:
    """
    Read a SnapGene (.dna) file and return:
      sequence  (str)  : DNA sequence in uppercase
      topology  (str)  : 'circular' | 'linear'
      name      (str)  : file name without extension
      features  (list) : [{name, type, start(0-based), end, strand(+/-), color}]
      primers   (list) : [{name, start(0-based), end, strand(+/-)}]
      notes     (dict) : free metadata from SnapGene
    """
    from xml.etree import ElementTree as ET

    with open(filepath, 'rb') as fh:
        data = fh.read()

    result = {
        'sequence': '', 'topology': 'linear',
        'name': Path(filepath).stem,
        'features': [], 'primers': [], 'notes': {},
    }

    for pkt_type, payload in _snapgene_read_packets(data):

        if pkt_type == 0 and len(payload) >= 2:
            flags = payload[0]
            result['topology'] = 'circular' if (flags & 0x01) else 'linear'
            result['sequence'] = payload[1:].decode('ascii', errors='replace').upper()

        elif pkt_type == 8:
            try:
                xml = ET.fromstring(payload.decode('utf-8', errors='replace'))
                for feat in xml.findall('.//Feature'):
                    name  = feat.get('name', 'unknown')
                    ftype = feat.get('type', 'misc_feature')
                    color = feat.get('color', '#aaaaaa')
                    direc = feat.get('directionality', '1')
                    strand = '+' if direc in ('1', 'forward') else '-'
                    for seg in feat.findall('Segment'):
                        rng = seg.get('range', '')
                        if '-' in rng:
                            s, e = rng.split('-', 1)
                            try:
                                result['features'].append({
                                    'name': name, 'type': ftype, 'color': color,
                                    'start': int(s) - 1,
                                    'end':   int(e),
                                    'strand': strand,
                                })
                            except ValueError:
                                pass
            except ET.ParseError:
                pass

        elif pkt_type == 6:
            try:
                xml = ET.fromstring(payload.decode('utf-8', errors='replace'))
                for pr in xml.findall('.//Primer'):
                    pname   = pr.get('name', 'primer')
                    pstrand = '+' if pr.get('templateStrand', 'sense') == 'sense' else '-'
                    for bind in pr.findall('BindingSite'):
                        rng = bind.get('location', '')
                        if '-' in rng:
                            s, e = rng.split('-', 1)
                            try:
                                result['primers'].append({
                                    'name': pname,
                                    'start': int(s) - 1, 'end': int(e),
                                    'strand': pstrand,
                                })
                            except ValueError:
                                pass
            except ET.ParseError:
                pass

        elif pkt_type == 10:
            try:
                xml = ET.fromstring(payload.decode('utf-8', errors='replace'))
                for child in xml:
                    result['notes'][child.tag] = (child.text or '').strip()
            except ET.ParseError:
                pass

    return result


def _sg_packet(pkt_type: int, payload: bytes) -> bytes:
    """Serialise a SnapGene TLV packet."""
    return bytes([pkt_type]) + _struct.pack('>I', len(payload)) + payload


def write_snapgene_dna(filepath: str, sequence: str, features: list,
                       primers: list = None, topology: str = 'linear',
                       name: str = '', notes: dict = None):
    """
    Write a SnapGene-compatible binary .dna file.

    Parameters
    ----------
    sequence  : str   — DNA in uppercase
    features  : list  — [{name, type, start(0-based), end, strand(+/-), color}]
    primers   : list  — [{name, start(0-based), end, strand(+/-)}]  (optional)
    topology  : 'linear' | 'circular'
    name      : str   — record name (written to Notes packet)
    notes     : dict  — additional metadata
    """
    from xml.etree.ElementTree import Element, SubElement, tostring as _tostr

    buf = _io.BytesIO()
    buf.write(b'\x09')   # magic cookie

    # ── Packet 0: sequence ──────────────────────────────────────
    flags = 0x01 if topology == 'circular' else 0x00
    seq_payload = bytes([flags]) + sequence.lower().encode('ascii')
    buf.write(_sg_packet(0, seq_payload))

    # ── Packet 8: features ──────────────────────────────────────
    if features:
        root_xml = Element('Features')
        for feat in features:
            direc = '1' if feat.get('strand', '+') == '+' else '2'
            f_el = SubElement(root_xml, 'Feature',
                              name=feat.get('name', 'feature'),
                              type=feat.get('type', 'misc_feature'),
                              directionality=direc,
                              color=feat.get('color', '#aaaaaa'),
                              swappedSegmentNumbering='0',
                              allowSegmentOverlaps='0')
            # SnapGene uses 1-based closed ranges
            s = feat['start'] + 1
            e = feat['end']
            SubElement(f_el, 'Segment', range=f"{s}-{e}",
                       name=feat.get('name', ''),
                       color=feat.get('color', '#aaaaaa'),
                       type='standard')
        xml_bytes = _tostr(root_xml, encoding='unicode').encode('utf-8')
        buf.write(_sg_packet(8, xml_bytes))

    # ── Packet 6: primers ───────────────────────────────────────
    primers = primers or []
    if primers:
        root_xml = Element('Primers')
        for pr in primers:
            tpl_strand = 'sense' if pr.get('strand', '+') == '+' else 'antisense'
            pr_el = SubElement(root_xml, 'Primer',
                               name=pr.get('name', 'primer'),
                               templateStrand=tpl_strand)
            s = pr['start'] + 1
            e = pr['end']
            SubElement(pr_el, 'BindingSite', location=f"{s}-{e}")
        xml_bytes = _tostr(root_xml, encoding='unicode').encode('utf-8')
        buf.write(_sg_packet(6, xml_bytes))

    # ── Packet 10: notes ────────────────────────────────────────
    notes = notes or {}
    if name:
        notes.setdefault('Description', name)
    notes.setdefault('CustomMapLabel', name or 'ORF Pipeline v20')
    notes.setdefault('ConfirmedExperimentally', '0')
    notes.setdefault('Created', datetime.now().strftime('%Y-%m-%d'))
    notes.setdefault('CreatedBy', 'ORF Secretion Pipeline v20')
    root_xml = Element('Notes')
    for k, v in notes.items():
        el = SubElement(root_xml, k)
        el.text = str(v)
    xml_bytes = _tostr(root_xml, encoding='unicode').encode('utf-8')
    buf.write(_sg_packet(10, xml_bytes))

    with open(filepath, 'wb') as fh:
        fh.write(buf.getvalue())
