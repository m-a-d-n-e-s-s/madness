#include <WorkflowBuilders.hpp>

#include <array>
#include <iostream>
#include <string>
#include <vector>

// Unit test for the madqc workflow dispatch surface (name -> WorkflowKind and
// the runnable-workflow listing). The v2 StateParallelPlanner half of this test
// was removed with the v2 engine (M1 decoupling Stage 3); state-parallel
// planning now lives in molresponse_v3's subworld layer (doc 32) and is
// covered by test_fd_subworld_fanout / test_node_subworlds.

int main() {
  using madness::workflow_builders::WorkflowKind;
  using madness::workflow_builders::workflow_kind_from_name;

  struct Case {
    const char *name;
    WorkflowKind kind;
  };

  constexpr std::array<Case, 8> expected = {{
      {"scf", WorkflowKind::Scf},
      {"nemo", WorkflowKind::Nemo},
      {"response", WorkflowKind::Response},
      {"mp2", WorkflowKind::Mp2Cc2},
      {"cc2", WorkflowKind::Mp2Cc2},
      {"cis", WorkflowKind::Cis},
      {"oep", WorkflowKind::Oep},
      {"optimize", WorkflowKind::Optimize},
  }};

  bool ok = true;
  for (const auto &entry : expected) {
    const auto actual = workflow_kind_from_name(entry.name);
    if (actual != entry.kind) {
      ok = false;
      std::cerr << "workflow_kind_from_name mismatch for '" << entry.name
                << "'\n";
    }
  }

  if (workflow_kind_from_name("unknown_workflow") != WorkflowKind::Unknown) {
    ok = false;
    std::cerr << "workflow_kind_from_name should return Unknown for invalid "
                 "workflow\n";
  }

  const std::string workflow_list =
      madness::workflow_builders::runnable_workflow_list();
  // `optimize` maps to a WorkflowKind (so --wf=optimize gets a migration message
  // rather than "unknown") but must NOT be advertised as a runnable workflow: it
  // is a flag on one, --optimize --wf=<scf|nemo>.
  if (workflow_list.find("optimize") != std::string::npos) {
    ok = false;
    std::cerr << "runnable_workflow_list must not advertise 'optimize'\n";
  }
  for (const char *name : madness::workflow_builders::runnable_workflows) {
    if (workflow_list.find(name) == std::string::npos) {
      ok = false;
      std::cerr << "runnable_workflow_list missing '" << name << "'\n";
    }
  }

  // --- StepContext: the nuclear-framework comparison the consumer-side guard
  // and the checkpoint guard both rely on. Molecule needs no World, so this
  // stays a plain unit test with no madness::initialize.
  {
    using madness::Molecule;
    using madness::same_nuclear_framework;

    auto make = [](double z1, double z2, int za = 8, bool third = false) {
      Molecule m;
      m.add_atom(0.0, 0.0, z1, static_cast<double>(za), za);
      m.add_atom(0.0, 0.0, z2, 1.0, 1);
      if (third)
        m.add_atom(0.0, 0.0, 5.0, 1.0, 1);
      return m;
    };

    struct Case {
      const char *what;
      Molecule a, b;
      bool expected;
    };
    const std::vector<Case> cases = {
        {"identical", make(0.0, 1.8), make(0.0, 1.8), true},
        // below Atom::operator=='s 1e-10 displacement threshold: JSON
        // round-trip / re-orientation noise must not count as a difference
        {"1e-12 displacement", make(0.0, 1.8), make(0.0, 1.8 + 1e-12), true},
        // a real geometry change (an optimizer step) must count
        {"1e-3 displacement", make(0.0, 1.8), make(0.0, 1.8 + 1e-3), false},
        {"different element", make(0.0, 1.8), make(0.0, 1.8, 7), false},
        {"extra atom", make(0.0, 1.8), make(0.0, 1.8, 8, true), false},
    };

    for (const auto &c : cases) {
      if (same_nuclear_framework(c.a, c.b) != c.expected) {
        ok = false;
        std::cerr << "same_nuclear_framework wrong for case '" << c.what
                  << "': expected " << c.expected << "\n";
      }
    }

    // an empty context provides nothing -- both reference slots unset, so a
    // consumer's guard is a no-op and a first/standalone step sees no artifacts
    const madness::StepContext empty;
    if (empty.reference || empty.nemo_reference || empty.molecule.has_value() ||
        !empty.archives.empty()) {
      ok = false;
      std::cerr << "default-constructed StepContext is not empty\n";
    }
  }

  return ok ? 0 : 1;
}
