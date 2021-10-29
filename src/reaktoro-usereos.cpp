#include "reaktoro-usereos.hpp"

#include <Reaktoro/Reaktoro.hpp>
#include <ThermoFun/ThermoFun.h>

#include <array>
#include <unordered_map>
#include <vector>

namespace {

const char *moduleName   = "GEMS Test module";
const int filenameLength = 256;
const int titleLength    = 80;

double fluidViscosity = 1e-3;
double solidViscosity = 1e7;

} // namespace

std::vector<double> linspace(double begin, double end, int count) {
  std::vector<double> result(count);
  for (int i = 0; i < count; ++i) {
    double x  = double(i) / (count - 1);
    result[i] = x * (end - begin) + begin;
  }
  return result;
}

class UserEosModule {
public:
  UserEosModule() : m_problem(m_system) {}

  int initialize(const std::string &configPath);
  void getParameters(char cmpNames[], double molWeights[], char phNames[], char auxNames[]);
  bool calculateEquilibrium(const double P, const double T, const double z[], int &nPhase, double props[],
                            int8_t phaseId[], double auxArray[], int8_t mode);

  int nComponents() const;
  int nMaxPhases() const;
  int nEndMembers() const;

private:
  void updateParams();
  void logFailedParams(const double P, const double T, const double z[]) const;

  Reaktoro::ChemicalEditor m_editor;
  Reaktoro::ChemicalSystem m_system;

  Reaktoro::EquilibriumProblem m_problem;
  Reaktoro::EquilibriumSolver m_solver;

  std::vector<int> m_fluidPhasesIds;
  std::vector<int> m_solidPhasesIds;

  std::vector<std::array<char, 3>> m_componentNames;
  std::vector<std::array<char, 3>> m_phasesNames;
  std::vector<std::array<char, 8>> m_speciesNames;

  std::vector<double> m_molarWeights;

  int m_numTransportPhases = 0;
  bool m_hasSolidPhases    = false;
};

int UserEosModule::initialize(const std::string &configPath) {
  std::cout << "Loading database... ";
  ThermoFun::Database db;
  try {
    db       = ThermoFun::Database(configPath);
    m_editor = Reaktoro::ChemicalEditor(db);
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }

  std::cout << "Done" << std::endl;

  std::vector<double> pressures    = linspace(1000, 2500, 51);
  std::vector<double> temperatures = linspace(200, 500, 51);

  m_editor.setPressures(pressures, "bar");
  m_editor.setTemperatures(temperatures, "degC");

  try {
    m_editor.addAqueousPhaseWithElements("H O Na Cl Cu S Si");
    m_editor.addMineralPhase("Halite Chalcocite Quartz");
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }

  m_system  = Reaktoro::ChemicalSystem(m_editor);
  m_problem = Reaktoro::EquilibriumProblem(m_system);
  m_solver  = Reaktoro::EquilibriumSolver(m_system);

  updateParams();

  return 0;
}

void UserEosModule::updateParams() {
  m_componentNames.clear();
  m_phasesNames.clear();
  m_speciesNames.clear();

  m_molarWeights.clear();

  m_fluidPhasesIds.clear();
  m_solidPhasesIds.clear();

  m_numTransportPhases = 0;
  m_hasSolidPhases     = false;

  const auto &phases = m_system.phases();
  for (int idx = 0; idx < phases.size(); ++idx) {
    const auto &phase = phases[idx];
    if (phase.isFluid()) {
      m_numTransportPhases++;
      m_fluidPhasesIds.push_back(idx);

      std::array<char, 3> name = {' ', ' ', ' '};
      std::copy_n(phase.name().data(), std::min(phase.name().length(), size_t(3)), name.begin());
      m_phasesNames.push_back(name);

      std::cout << "  New fluid phase: " << phase.name() << std::endl;
    } else if (phase.isSolid()) {
      if (!m_hasSolidPhases) {
        m_numTransportPhases++;
        m_hasSolidPhases = true;
      }
      m_solidPhasesIds.push_back(idx);
    }
  }

  if (m_hasSolidPhases) {
    m_phasesNames.push_back({'S', 'O', 'L'});
  }

  std::cout << "Elements: ";
  for (const auto &elem : m_system.elements()) {
    std::array<char, 3> name = {' ', ' ', ' '};
    std::copy_n(elem.name().data(), std::min(elem.name().length(), size_t(3)), name.begin());
    m_componentNames.push_back(name);
    m_molarWeights.push_back(elem.molarMass());

    std::cout << elem.name() << ' ';
  }
  std::cout << std::endl;

  std::unordered_map<std::string, int> speciesMap;
  std::string charsToRemove = "+-#@()";

  std::cout << "Species: ";
  for (const auto &specie : m_system.species()) {
    std::array<char, 8> name = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};
    std::string aux("#");
    std::string esc = specie.name().substr(0, 7);

    if (speciesMap.find(esc) == speciesMap.end()) {
      speciesMap[esc] = 0;
    } else {
      speciesMap[esc]++;
    }
    auto prefix = std::to_string(speciesMap[esc]);
    aux.append(esc);
    std::copy_n(aux.data(), std::min(aux.length(), size_t(8)), name.begin());
    m_speciesNames.push_back(name);

    std::cout << specie.name() << " (" << aux << ") ";
  }
  std::cout << std::endl;
}

void UserEosModule::getParameters(char cmpNames[], double molWeights[], char phNames[], char auxNames[]) {
  memcpy(cmpNames, m_componentNames.data(), 3 * size_t(nComponents()));
  memcpy(phNames, m_phasesNames.data(), 3 * size_t(nMaxPhases()));
  if (auxNames != nullptr) {
    memcpy(auxNames, m_speciesNames.data(), 8 * size_t(nEndMembers()));
  }
  std::copy(m_molarWeights.begin(), m_molarWeights.end(), molWeights);
}

bool UserEosModule::calculateEquilibrium(const double P, const double T, const double z[], int &nPhase, double props[],
                                         int8_t phaseId[], double auxArray[], int8_t mode) {
  m_problem.setPressure(P);
  m_problem.setTemperature(T);

  double molesTotal = 0.0;
  for (int ic = 0; ic < nComponents(); ++ic) {
    m_problem.setElementAmount(ic, z[ic] / m_molarWeights[ic]);
  }

  Reaktoro::ChemicalState state(m_system);
  Reaktoro::EquilibriumResult result;

  try {
    result = m_solver.solve(state, m_problem);
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "System params:" << std::endl;
    logFailedParams(P, T, z);
  }

  if (!result.optimum.succeeded) {
    std::cerr << "Error: the problem did not converge, system params:" << std::endl;
    logFailedParams(P, T, z);
  }

  size_t offset = size_t(6 + nComponents());
  std::fill(props, props + nMaxPhases() * offset, 0.0);

  auto properties = state.properties();

  auto phaseMasses     = properties.phaseMasses();
  auto phaseVolumes    = properties.phaseVolumes();
  auto phaseDensities  = properties.phaseDensities();
  auto phaseEnthalpies = properties.phaseSpecificEnthalpies();

  auto totalVolume = double(properties.volume());

  for (int ip = 0; ip < m_fluidPhasesIds.size(); ++ip) {
    int iphase = m_fluidPhasesIds[ip];

    props[ip * offset + 0] = double(phaseDensities[iphase]);
    props[ip * offset + 1] = double(phaseEnthalpies[iphase]);
    props[ip * offset + 2] = fluidViscosity;
    props[ip * offset + 3] = double(phaseVolumes[iphase]) / totalVolume;

    double sumC = 0.0;

    for (int ic = 0; ic < nComponents(); ++ic) {
      double x                    = state.elementAmountInPhase(ic, iphase);
      double phaseMass            = double(phaseMasses[iphase]);
      props[ip * offset + 6 + ic] = x * m_molarWeights[ic] / phaseMass;
      sumC += props[ip * offset + 6 + ic];
    }

    if (std::abs(sumC - 1.0) > 1e-15) {
      std::cerr << "Error: fluid concentrations do not sum to 1!!!\n";
    }

    phaseId[ip] = int8_t(ip + 1);
  }

  if (m_hasSolidPhases) {
    int ip = m_fluidPhasesIds.size();

    double solidVolume = double(properties.solidVolume());

    double solidMass     = 0.0;
    double solidEnthalpy = 0.0;

    for (int is = 0; is < m_solidPhasesIds.size(); ++is) {
      int iphase       = m_solidPhasesIds[is];
      double phaseMass = double(phaseMasses[iphase]);
      solidMass += phaseMass;
      solidEnthalpy += double(phaseEnthalpies[iphase]) * phaseMass;

      for (int ic = 0; ic < nComponents(); ++ic) {
        double x = state.elementAmountInPhase(ic, iphase);
        props[ip * offset + 6 + ic] += x * m_molarWeights[ic];
      }
    }

    // Specific enthalpy
    solidEnthalpy /= solidMass;

    props[ip * offset + 0] = solidMass / solidVolume;
    props[ip * offset + 1] = solidEnthalpy;
    props[ip * offset + 2] = solidViscosity;
    props[ip * offset + 3] = solidVolume / totalVolume;

    double sumC = 0.0;
    // Normalize the concentrations
    for (int ic = 0; ic < nComponents(); ++ic) {
      props[ip * offset + 6 + ic] /= solidMass;
      sumC += props[ip * offset + 6 + ic];
    }
    if (std::abs(sumC - 1.0) > 1e-15) {
      std::cerr << "Error: solid concentrations do not sum to 1!!!\n";
    }

    phaseId[ip] = int8_t(ip + 1);
  }

  nPhase = nMaxPhases();

  for (int iem = 0; iem < nEndMembers(); ++iem) {
    auxArray[iem] = state.speciesAmount(iem);
  }

  return true;
}

void UserEosModule::logFailedParams(const double P, const double T, const double z[]) const {
  std::cerr << "P: " << P / 1e5 << " bar\n";
  std::cerr << "T: " << T << " K\n";
  std::cerr << "Elemental concentrations:\n";
  for (int ic = 0; ic < nComponents(); ++ic) {
    std::cerr << m_system.element(ic).name() << ": " << z[ic] << '\n';
  }
}

int UserEosModule::nComponents() const { return m_system.numElements() - 1; }

int UserEosModule::nMaxPhases() const { return m_numTransportPhases; }

int UserEosModule::nEndMembers() const { return m_system.numSpecies(); }

static UserEosModule userModule;

void ReadConfigurationFile(char *filename, int *ierr, char *title) {
  std::string cFilename(filename, filename + filenameLength);
  size_t spacePos = cFilename.find(' ');
  if (spacePos != std::string::npos) {
    cFilename[cFilename.find(' ')] = '\0';
  }
  memset(title, ' ', size_t(titleLength));
  memcpy(title, moduleName, strlen(moduleName));

  try {
    *ierr = userModule.initialize(cFilename);
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}

void GetDimensions(int *nComponents, int *nPhaseMax, int *nAux) {
  *nComponents = userModule.nComponents();
  *nPhaseMax   = userModule.nMaxPhases();
  *nAux        = userModule.nEndMembers();

  std::cout << "nMaxPhases: " << userModule.nMaxPhases() << std::endl;
}

void GetGlobalParameters(char *cmpNames, double *molWeights, char *phNames, char *auxNames, char *auxUnits,
                         int8_t *opt) {
  try {
    userModule.getParameters(cmpNames, molWeights, phNames, auxNames);
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  char unit[] = "NODIM   ";
  for (int idx = 0; idx < userModule.nEndMembers(); ++idx) {
    memcpy(auxUnits + size_t(idx) * 8, unit, 8);
  }
  opt[0] = 0;
  opt[1] = 0;
}

void PhaseEquilibrium(double *pres, double *temp, double *z, int *nPhase, double *props, int8_t *phaseId,
                      double *auxArray, int8_t *mode) {
  userModule.calculateEquilibrium(*pres, *temp, z, *nPhase, props, phaseId, auxArray, *mode);
}
