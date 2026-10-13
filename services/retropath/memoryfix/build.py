"""Compile the pinned cleanup repair and allocator agent inside the Docker build."""
from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile

VERSION = "4.9.1-glade-memoryfix2"
KNIME = Path("/opt/knime/knime_4.7.0")
PLUGIN = KNIME / "p2/plugins/org.rdkit.knime.nodes_4.9.1.v202312060020.jar"
ORIGINAL_SHA = "bc42ad1cdda99718c442824414752a09389c617c8934b17c27680ee664e50e78"
EXPECTED_CLASS_SHA = "5fae50ef34b92bf68e0d3fb93740c91ae3a900f67c6f7a30196794bf1afc29de"
CLASS = "org/rdkit/knime/nodes/onecomponentreaction2/AbstractRDKitReactionNodeModel.class"
JAVA = KNIME / "plugins/org.knime.binary.jre.linux.x86_64_17.0.5.20221116/jre/bin/java"
COMPILER = KNIME / "plugins/org.eclipse.jdt.core_3.30.0.v20220524-1333.jar"
ASM = KNIME / "plugins/org.objectweb.asm_9.3.0.v20220409-0157.jar"
JNA = KNIME / "plugins/com.sun.jna_5.8.0.v20210503-0343.jar"
BASE = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(arguments: list[str], log: Path) -> None:
    environment = os.environ.copy()
    environment["_JAVA_OPTIONS"] = "-Xmx128m"
    environment["JAVA_TOOL_OPTIONS"] = "-XX:ActiveProcessorCount=1"
    result = subprocess.run([str(JAVA), *arguments], env=environment, capture_output=True, text=True, timeout=120)
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(log.read_text())


def main() -> None:
    if sha(PLUGIN) != ORIGINAL_SHA:
        raise RuntimeError("Unsupported original RDKit plugin; refusing to modify it")
    for dependency in (JAVA, COMPILER, ASM, JNA):
        if not dependency.is_file():
            raise RuntimeError("Pinned build dependency missing: " + str(dependency))
    with tempfile.TemporaryDirectory(prefix="glade-memoryfix-") as temporary:
        work = Path(temporary)
        with zipfile.ZipFile(PLUGIN) as source:
            (work / "original.class").write_bytes(source.read(CLASS))
        run(["-cp", str(COMPILER), "org.eclipse.jdt.internal.compiler.batch.Main", "-proc:none", "-17",
            "-cp", str(ASM) + ":" + str(JNA), "-d", str(work), str(BASE / "RegisterReactionVector.java"),
            str(BASE / "NativeMemoryAgent.java")], BASE / "compile.log")
        run(["-cp", str(work) + ":" + str(ASM), "RegisterReactionVector", str(work / "original.class"),
            str(work / "patched.class")], BASE / "patch.log")
        if sha(work / "patched.class") != EXPECTED_CLASS_SHA:
            raise RuntimeError("Generated node class differs from the chemically validated repair")
        removed = []
        output = work / "patched.jar"
        with zipfile.ZipFile(PLUGIN) as source, zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as target:
            for info in source.infolist():
                if info.filename.startswith("META-INF/") and info.filename.upper().endswith((".SF", ".RSA", ".DSA", ".EC")):
                    removed.append(info.filename)
                    continue
                data = (work / "patched.class").read_bytes() if info.filename == CLASS else source.read(info.filename)
                target.writestr(copy.copy(info), data)
        with zipfile.ZipFile(PLUGIN) as original, zipfile.ZipFile(output) as patched:
            expected = set(original.namelist()) - set(removed)
            if set(patched.namelist()) != expected or len(expected) - 1 != 734:
                raise RuntimeError("Unexpected plugin member inventory")
            if any(original.read(name) != patched.read(name) for name in expected if name != CLASS):
                raise RuntimeError("A plugin member other than the intended node class changed")
        manifest = "Manifest-Version: 1.0\r\nPremain-Class: NativeMemoryAgent\r\nClass-Path: file:" + str(JNA) + "\r\n\r\n"
        lines = []
        for line in manifest.split("\r\n"):
            while len(line.encode()) > 70:
                lines.append(line[:70]); line = " " + line[70:]
            lines.append(line)
        classes = ("NativeMemoryAgent.class", "NativeMemoryAgent$MallInfo.class", "NativeMemoryAgent$LibC.class")
        with zipfile.ZipFile(BASE / "native-memory-agent.jar", "w", zipfile.ZIP_DEFLATED) as target:
            for name, data in [("META-INF/MANIFEST.MF", "\r\n".join(lines).encode()),
                    *((name, (work / name).read_bytes()) for name in classes)]:
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                target.writestr(info, data)
        output.replace(PLUGIN)
        record = {"version": VERSION, "original_jar_sha256": ORIGINAL_SHA, "original_plugin_path": str(PLUGIN),
            "modified_member": CLASS, "patched_class_sha256": EXPECTED_CLASS_SHA, "removed_signatures": removed,
            "unchanged_members": 734, "build_method": "Compile Java repair/agent from source during Docker build; modify one pinned installed node class",
            "source_sha256": {name: sha(BASE / name) for name in ("RegisterReactionVector.java", "NativeMemoryAgent.java", "rdkit-vector-cleanup.patch", "build.py")},
            "agent_class_sha256": {name: sha(work / name) for name in classes},
            "artifacts_sha256": {"rdkit-nodes-memoryfix.jar": sha(PLUGIN), "native-memory-agent.jar": sha(BASE / "native-memory-agent.jar")},
            "rules_sha256": "e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f",
            "chemistry": "No changes to rules, explicit hydrogen, enumeration limits, sanitization, deduplication or scoring"}
        (BASE / "patch_manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
