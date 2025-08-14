# kegg_module.py
"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""

from Bio.KEGG import REST
from Bio.KEGG.Enzyme import parse as kegg_enzyme_parse
import io
import requests

KEGG_BASE = "https://rest.kegg.jp"

def _ensure_ec_id(ec_number: str) -> str:
    ec_raw = (ec_number or "").strip()
    return ec_raw if ec_raw.lower().startswith("ec:") else f"ec:{ec_raw}"

def _pathways_from_enzyme_record(record):
    pathways = []
    for p in getattr(record, "pathway", []) or []:
        if isinstance(p, tuple) and len(p) >= 3:
            db, pid, pname = p[0], p[1], p[2]
            pathways.append({"db": str(db), "id": str(pid), "name": str(pname)})
    return pathways

def _link_ec_to_pathways(ec_number: str):
    try:
        url = f"{KEGG_BASE}/link/pathway/ec:{ec_number.replace('ec:', '')}"
        r = requests.get(url, timeout=20)
        if r.status_code != 200 or not r.text.strip():
            return []
        ids = []
        for line in r.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2 and parts[1].startswith("path:"):
                ids.append(parts[1])
        ids = sorted(set(ids))
        if not ids:
            return []
        list_q = "+".join(ids)
        lr = requests.get(f"{KEGG_BASE}/list/{list_q}", timeout=20)
        if lr.status_code != 200 or not lr.text.strip():
            return [{"db": "path", "id": pid, "name": ""} for pid in ids]
        name_map = {}
        for line in lr.text.strip().split("\n"):
            if "\t" in line:
                pid, desc = line.split("\t", 1)
                name_map[pid] = desc
        return [{"db": "path", "id": pid, "name": name_map.get(pid, "")} for pid in ids]
    except Exception:
        return []

def fetch_kegg_enzyme_info(ec_number):
    try:
        kegg_id = _ensure_ec_id(ec_number)
        response = REST.kegg_get(kegg_id)
        data = response.read()
        if not data or not data.strip():
            raise ValueError("Empty KEGG response")
        record = list(kegg_enzyme_parse(io.StringIO(data)))[0]
        pathways = _pathways_from_enzyme_record(record)
        if not pathways:
            pathways = _link_ec_to_pathways(record.entry or ec_number)
        return {
            "status": "success",
            "type": "enzyme",
            "ec_number": getattr(record, "entry", ec_number),
            "name": getattr(record, "name", []),
            "class": getattr(record, "classname", []),
            "pathways": pathways,
            "cofactors": getattr(record, "cofactor", []),
            "effectors": getattr(record, "effector", []),
        }
    except Exception as e:
        return {"status": "error", "error_message": str(e)}

def search_kegg_pathway(query):
    try:
        url = f"{KEGG_BASE}/find/pathway/{query}"
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            return {"status": "error", "error_message": f"KEGG search failed ({r.status_code})"}
        results = []
        for line in r.text.strip().split("\n"):
            if "\t" not in line:
                continue
            kegg_id, description = line.split("\t", 1)
            results.append({"id": kegg_id, "description": description})
        return {"status": "success", "results": results}
    except Exception as e:
        return {"status": "error", "error_message": str(e)}

def fetch_kegg_details(kegg_id):
    try:
        url = f"{KEGG_BASE}/get/{kegg_id}"
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            return {"status": "error", "error_message": f"KEGG details fetch failed ({r.status_code})"}
        return {"status": "success", "details": r.text}
    except Exception as e:
        return {"status": "error", "error_message": str(e)}

def search_kegg(query):
    """Alias for backward compatibility with old imports."""
    return search_kegg_pathway(query)
