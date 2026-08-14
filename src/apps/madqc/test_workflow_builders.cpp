#include <WorkflowBuilders.hpp>

#include <array>
#include <iostream>
#include <string>

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

  return ok ? 0 : 1;
}
