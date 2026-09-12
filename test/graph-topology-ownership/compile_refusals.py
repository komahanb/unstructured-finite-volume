"""Reject direct mutation and pointer access through the public topology API."""

from pathlib import Path
import subprocess
import sys
import tempfile


def source(statement, declarations="", procedures=""):
    return f"""program topology_contract
  use view_directed_stored, only: stored_directed_graph
  use relation_binary, only: csr_relation, integer_fibre
  implicit none
  type(stored_directed_graph), target :: g
  type(csr_relation), target :: r
  type(integer_fibre) :: f
  {declarations}
  {statement}
  {procedures}
end program topology_contract
"""


def main():
    compiler = sys.argv[1:]
    if not compiler:
        raise SystemExit("compiler and module search flags are required")
    refusals = {}
    for component in ("number", "nv", "ne", "reversed", "tail", "head", "num_without_head",
                      "xinc", "einc", "xadj", "vadj", "xout", "eout", "xin", "ein",
                      "vtag", "etag", "whole_rel", "vset", "eset"):
        refusals[component] = source(f"g % {component} = g % {component}")
    refusals["fibre_storage"] = source("f % entries(1) = 1")
    refusals["fibre_pointer"] = source("p => r % image_view(1)", "integer, pointer :: p(:)")
    refusals["fibre_member"] = source("f % member(1) = 1")
    refusals["fibre_mutable_reader"] = source("call f % read(reader)", procedures="""contains
  subroutine reader(members)
    integer, intent(inout) :: members(:)
    members = 1
  end subroutine reader""")
    refusals["graph_mutable_reader"] = source("call g % read_incoming(reader)", procedures="""contains
  subroutine reader(offsets, indices, sources)
    integer, intent(inout) :: offsets(:), indices(:), sources(:)
    sources = 1
  end subroutine reader""")
    refusals["fibre_target_reader"] = source("call f % read(reader)", procedures="""contains
  subroutine reader(members)
    integer, target, intent(in) :: members(:)
  end subroutine reader""")
    refusals["graph_target_reader"] = source("call g % read_incoming(reader)", procedures="""contains
  subroutine reader(offsets, indices, sources)
    integer, target, intent(in) :: offsets(:), indices(:), sources(:)
  end subroutine reader""")
    refusals["context_mutable_members"] = source("call f % read(reader, accumulator)",
                                                "integer :: accumulator", """contains
  subroutine reader(members, context)
    integer, intent(inout) :: members(:)
    class(*), intent(inout) :: context
    members = 1
  end subroutine reader""")
    refusals["context_target_members"] = source("call f % read(reader, accumulator)",
                                               "integer :: accumulator", """contains
  subroutine reader(members, context)
    integer, target, intent(in) :: members(:)
    class(*), intent(inout) :: context
  end subroutine reader""")
    refusals["context_wrong_intent"] = source("call f % read(reader, accumulator)",
                                             "integer :: accumulator", """contains
  subroutine reader(members, context)
    integer, intent(in) :: members(:)
    class(*), intent(in) :: context
  end subroutine reader""")
    refusals["context_target_argument"] = source("call f % read(reader, accumulator)",
                                                "integer :: accumulator", """contains
  subroutine reader(members, context)
    integer, intent(in) :: members(:)
    class(*), target, intent(inout) :: context
  end subroutine reader""")
    with tempfile.TemporaryDirectory(prefix="ufvm-topology-refusals-") as directory:
        path = Path(directory) / "contract.f90"
        path.write_text(source("g = stored_directed_graph(1, [1], [1])\n  print *, g % num_edges()\n"
                               "call f % read(reader, accumulator)", "integer :: accumulator", """contains
  subroutine reader(members, context)
    integer, intent(in) :: members(:)
    class(*), intent(inout) :: context
  end subroutine reader"""))
        control = subprocess.run(compiler + ["-fsyntax-only", str(path)], capture_output=True, text=True)
        if control.returncode != 0:
            raise SystemExit("valid public reader failed to compile:\n" + control.stderr)
        for name, program in refusals.items():
            path.write_text(program)
            result = subprocess.run(compiler + ["-fsyntax-only", str(path)], capture_output=True, text=True)
            if result.returncode == 0:
                raise SystemExit(f"FAIL : external mutation admitted: {name}")
            if name in ("fibre_storage", *tuple(refusals)[:20]) and "private" not in result.stderr.lower():
                raise SystemExit(f"FAIL : {name} refused for an unrelated reason:\n{result.stderr}")
            print(f"PASS : compiler refuses {name}")


if __name__ == "__main__":
    main()
