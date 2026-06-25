process SOFTWARE_VERSIONS_HTML {
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}"
    
    input:
        path(versions_files, stageAs: 'versions??/*')

    output:
        path "software_versions.html"
        path "software_versions.tsv"
        path "versions.yml"

    shell:
    '''
    python3 <<'PY'
    import re
    import html
    from pathlib import Path

    rows = []
    seen = set()

    for f in sorted(Path(".").glob("versions*/versions.yml")):
        text = f.read_text()
        current_process = None

        for line in text.splitlines():
            if not line.strip():
                continue

            if re.match(r'^".*":\\s*$', line) or re.match(r'^[^\\s].*:\\s*$', line):
                current_process = line.strip().rstrip(':').strip('"').strip("'")
                continue

            m = re.match(r'^\\s+([^:]+):\\s*(.*)$', line)
            if m and current_process:
                tool = m.group(1).strip()
                version = m.group(2).strip().strip('"').strip("'")

                key = (current_process, tool, version)
                if key not in seen:
                    seen.add(key)
                    rows.append((current_process, tool, version))

    rows.sort(key=lambda x: (x[0].lower(), x[1].lower(), x[2].lower()))

    with open("software_versions.tsv", "w") as out:
        out.write("process\\ttool\\tversion\\n")
        for process, tool, version in rows:
            out.write(f"{process}\\t{tool}\\t{version}\\n")

    with open("software_versions.html", "w") as out:
        out.write("<html><head><title>Software Versions</title></head><body>")
        out.write("<h1>Software Versions</h1>")
        out.write("<table border='1' cellpadding='5' cellspacing='0'>")
        out.write("<tr><th>Process</th><th>Tool</th><th>Version</th></tr>")

        for process, tool, version in rows:
            out.write(
                f"<tr><td>{html.escape(process)}</td>"
                f"<td>{html.escape(tool)}</td>"
                f"<td>{html.escape(version)}</td></tr>"
            )

        out.write("</table></body></html>")
    PY

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      python: $(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    '''

}