#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TRandom.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "nlohmann/json.hpp"

namespace fs = std::filesystem;

namespace {

using Garfield::AvalancheMicroscopic;
using Garfield::ComponentAnalyticField;
using Garfield::MediumMagboltz;
using Garfield::Sensor;
using json = nlohmann::json;

constexpr double kNmToCm = 1.e-7;
constexpr double kKvPerCmToVPerCm = 1.e3;
constexpr double kDefaultTemperatureK = 293.15;
constexpr double kDefaultInitialElectronEnergyEv = 0.1;
constexpr double kDefaultMaxElectronEnergyEv = 200.;

struct CliOptions {
  fs::path configPath{"config/default_parallel_plate.json"};
  fs::path outDir{"results"};
  std::string scanSpec{"all"};
};

struct BaselineConfig {
  double distanceNm = 300.;
  double pressureTorr = 760.;
  double fieldKvPerCm = 200.;
};

struct ScanGridConfig {
  std::vector<double> distanceNm{200., 250., 300., 350.};
  std::vector<double> pressureTorr{200., 400., 600., 760.};
  std::vector<double> fieldKvPerCm{50., 100., 200., 400.};
};

struct GasConfig {
  bool enablePenning = true;
  bool enableThermalMotion = true;
};

struct SimulationConfig {
  std::size_t numPrimaries = 500;
  double initialElectronEnergyEv = kDefaultInitialElectronEnergyEv;
  double maxElectronEnergyEv = kDefaultMaxElectronEnergyEv;
  std::size_t avalancheSizeLimit = 100000;
  double temperatureK = kDefaultTemperatureK;
};

struct HistogramConfig {
  int energyBins = 200;
  int multiplicityBinsMin = 20;
};

struct Config {
  BaselineConfig baseline;
  ScanGridConfig scan;
  GasConfig gas;
  SimulationConfig simulation;
  HistogramConfig histogram;
};

enum class ScanKind { Distance, Pressure, Field };

struct ScanPoint {
  ScanKind kind;
  std::string scanKey;
  std::string parameterKey;
  std::string parameterUnit;
  std::string xAxisTitle;
  double parameterValue = 0.;
  double distanceNm = 0.;
  double pressureTorr = 0.;
  double fieldKvPerCm = 0.;
  std::string tag;
  std::string label;
};

struct PointSummary {
  std::string scanKey;
  std::string parameterKey;
  std::string parameterUnit;
  double parameterValue = 0.;
  std::size_t numPrimaries = 0;
  std::size_t totalArrivals = 0;
  std::size_t primariesWithArrival = 0;
  std::size_t sizeLimitHitPrimaries = 0;
  std::size_t transportFailures = 0;
  double arrivalFraction = 0.;
  double meanEnergyEv = 0.;
  double rmsEnergyEv = 0.;
  double semEnergyEv = 0.;
  double meanMultiplicity = 0.;
  double rmsMultiplicity = 0.;
  double semMultiplicity = 0.;
  std::string xAxisTitle;
};

struct PointData {
  std::vector<double> arrivalEnergiesEv;
  std::vector<double> arrivalMultiplicities;
  PointSummary summary;
};

struct EndpointStatusBreakdown {
  std::size_t hitArrivalPlane = 0;
  std::size_t hitSourcePlane = 0;
  std::size_t hitOtherPlane = 0;
  std::size_t attached = 0;
  std::size_t belowTransportCut = 0;
  std::size_t leftDriftMediumArrivalBoundary = 0;
  std::size_t leftDriftMediumSourceBoundary = 0;
  std::size_t leftDriftMediumOther = 0;
  std::size_t leftDriftArea = 0;
  std::size_t outsideTimeWindow = 0;
  std::size_t calculationAbandoned = 0;
  std::size_t other = 0;
};

std::string JoinPath(const std::initializer_list<std::string_view>& parts,
                     const std::string_view separator) {
  std::ostringstream stream;
  bool first = true;
  for (const auto part : parts) {
    if (!first) stream << separator;
    first = false;
    stream << part;
  }
  return stream.str();
}

std::string ToLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

std::string FormatNumber(const double value, const int precision = 6) {
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(precision) << value;
  std::string text = stream.str();
  while (!text.empty() && text.back() == '0') text.pop_back();
  if (!text.empty() && text.back() == '.') text.pop_back();
  if (text.empty()) return "0";
  return text;
}

std::string FileSafeNumber(const double value) {
  std::string text = FormatNumber(value);
  std::replace(text.begin(), text.end(), '.', 'p');
  std::replace(text.begin(), text.end(), '-', 'm');
  return text;
}

std::string ShellQuote(const std::string_view value) {
  std::string quoted;
  quoted.reserve(value.size() + 2);
  quoted.push_back('"');
  for (const char c : value) {
    if (c == '"' || c == '\\') quoted.push_back('\\');
    quoted.push_back(c);
  }
  quoted.push_back('"');
  return quoted;
}

std::string BuildCommandLine(const int argc, char* argv[]) {
  std::ostringstream stream;
  for (int i = 0; i < argc; ++i) {
    if (i > 0) stream << ' ';
    stream << ShellQuote(argv[i]);
  }
  return stream.str();
}

[[noreturn]] void ThrowJsonTypeError(
    const std::initializer_list<std::string_view>& path,
    const std::string_view expectation) {
  throw std::runtime_error("Invalid JSON at '" + JoinPath(path, ".") +
                           "': expected " + std::string(expectation) + ".");
}

const json* FindObjectMember(const json& object, const std::string_view key,
                             const std::initializer_list<std::string_view>& path) {
  if (!object.is_object()) ThrowJsonTypeError(path, "an object");
  const auto iterator = object.find(std::string(key));
  if (iterator == object.end()) return nullptr;
  return &(*iterator);
}

const json* FindObjectSection(const json& object, const std::string_view key) {
  const auto* section = FindObjectMember(object, key, {key});
  if (!section) return nullptr;
  if (!section->is_object()) ThrowJsonTypeError({key}, "an object");
  return section;
}

double ReadOptionalDouble(const json& object, const std::string_view section,
                          const std::string_view key, const double fallback) {
  const auto* value = FindObjectMember(object, key, {section, key});
  if (!value) return fallback;
  if (!value->is_number()) ThrowJsonTypeError({section, key}, "a number");
  return value->get<double>();
}

std::size_t ReadOptionalSizeT(const json& object, const std::string_view section,
                              const std::string_view key,
                              const std::size_t fallback) {
  const auto* value = FindObjectMember(object, key, {section, key});
  if (!value) return fallback;
  if (!value->is_number_integer() && !value->is_number_unsigned()) {
    ThrowJsonTypeError({section, key}, "a non-negative integer");
  }
  const auto parsed = value->get<long long>();
  if (parsed < 0) {
    throw std::runtime_error("Invalid JSON at '" + JoinPath({section, key}, ".") +
                             "': expected a non-negative integer.");
  }
  return static_cast<std::size_t>(parsed);
}

int ReadOptionalInt(const json& object, const std::string_view section,
                    const std::string_view key, const int fallback) {
  const auto* value = FindObjectMember(object, key, {section, key});
  if (!value) return fallback;
  if (!value->is_number_integer()) {
    ThrowJsonTypeError({section, key}, "an integer");
  }
  return value->get<int>();
}

bool ReadOptionalBool(const json& object, const std::string_view section,
                      const std::string_view key, const bool fallback) {
  const auto* value = FindObjectMember(object, key, {section, key});
  if (!value) return fallback;
  if (!value->is_boolean()) ThrowJsonTypeError({section, key}, "a boolean");
  return value->get<bool>();
}

std::vector<double> ReadOptionalDoubleArray(const json& object,
                                            const std::string_view section,
                                            const std::string_view key,
                                            const std::vector<double>& fallback) {
  const auto* value = FindObjectMember(object, key, {section, key});
  if (!value) return fallback;
  if (!value->is_array()) ThrowJsonTypeError({section, key}, "an array of numbers");
  std::vector<double> values;
  values.reserve(value->size());
  for (std::size_t i = 0; i < value->size(); ++i) {
    const auto& item = value->at(i);
    if (!item.is_number()) {
      throw std::runtime_error("Invalid JSON at '" +
                               JoinPath({section, key, std::to_string(i)}, ".") +
                               "': expected a number.");
    }
    values.push_back(item.get<double>());
  }
  return values;
}

json ReadJsonFile(const fs::path& path) {
  std::ifstream stream(path);
  if (!stream) {
    throw std::runtime_error("Failed to open JSON file: " + path.string());
  }
  try {
    return json::parse(stream);
  } catch (const json::parse_error& error) {
    throw std::runtime_error("Failed to parse JSON file '" + path.string() +
                             "': " + error.what());
  }
}

std::string JoinPathTokenValues(const std::vector<double>& values,
                                const std::string_view unit) {
  std::ostringstream stream;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i > 0) stream << '-';
    stream << FileSafeNumber(values[i]);
  }
  stream << unit;
  return stream.str();
}

template <typename T>
double Mean(const std::vector<T>& values) {
  if (values.empty()) return 0.;
  const double sum = std::accumulate(values.begin(), values.end(), 0.0);
  return sum / static_cast<double>(values.size());
}

template <typename T>
double Rms(const std::vector<T>& values, const double mean) {
  if (values.empty()) return 0.;
  double variance = 0.;
  for (const auto value : values) {
    const double delta = static_cast<double>(value) - mean;
    variance += delta * delta;
  }
  variance /= static_cast<double>(values.size());
  return std::sqrt(variance);
}

double Sem(const double rms, const std::size_t n) {
  if (n == 0) return 0.;
  return rms / std::sqrt(static_cast<double>(n));
}

double NmToCm(const double valueNm) { return valueNm * kNmToCm; }
double KvPerCmToVPerCm(const double valueKvPerCm) {
  return valueKvPerCm * kKvPerCmToVPerCm;
}

void EnsureDirectory(const fs::path& path) { fs::create_directories(path); }

[[noreturn]] void PrintUsageAndExit(const char* program, const int exitCode) {
  std::ostream& out = exitCode == 0 ? std::cout : std::cerr;
  out << "Usage: " << program
      << " [--config <path>] [--scan all|d|p|field] [--out <directory>]\n";
  std::exit(exitCode);
}

CliOptions ParseCommandLine(const int argc, char* argv[]) {
  CliOptions options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--config") {
      if (i + 1 >= argc) PrintUsageAndExit(argv[0], 1);
      options.configPath = argv[++i];
    } else if (arg == "--scan") {
      if (i + 1 >= argc) PrintUsageAndExit(argv[0], 1);
      options.scanSpec = argv[++i];
    } else if (arg == "--out") {
      if (i + 1 >= argc) PrintUsageAndExit(argv[0], 1);
      options.outDir = argv[++i];
    } else if (arg == "--help" || arg == "-h") {
      PrintUsageAndExit(argv[0], 0);
    } else {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }
  return options;
}

Config LoadConfig(const fs::path& configPath) {
  if (!fs::exists(configPath)) {
    throw std::runtime_error("Configuration file not found: " + configPath.string());
  }

  const json root = ReadJsonFile(configPath);
  if (!root.is_object()) ThrowJsonTypeError({"<root>"}, "a JSON object");

  Config config;
  if (const auto* baseline = FindObjectSection(root, "baseline")) {
    config.baseline.distanceNm =
        ReadOptionalDouble(*baseline, "baseline", "distance_nm",
                           config.baseline.distanceNm);
    config.baseline.pressureTorr =
        ReadOptionalDouble(*baseline, "baseline", "pressure_torr",
                           config.baseline.pressureTorr);
    config.baseline.fieldKvPerCm =
        ReadOptionalDouble(*baseline, "baseline", "field_kv_per_cm",
                           config.baseline.fieldKvPerCm);
  }

  if (const auto* scan = FindObjectSection(root, "scan")) {
    config.scan.distanceNm =
        ReadOptionalDoubleArray(*scan, "scan", "distance_nm",
                                config.scan.distanceNm);
    config.scan.pressureTorr =
        ReadOptionalDoubleArray(*scan, "scan", "pressure_torr",
                                config.scan.pressureTorr);
    config.scan.fieldKvPerCm =
        ReadOptionalDoubleArray(*scan, "scan", "field_kv_per_cm",
                                config.scan.fieldKvPerCm);
  }

  if (const auto* gas = FindObjectSection(root, "gas")) {
    config.gas.enablePenning =
        ReadOptionalBool(*gas, "gas", "enable_penning",
                         config.gas.enablePenning);
    config.gas.enableThermalMotion =
        ReadOptionalBool(*gas, "gas", "enable_thermal_motion",
                         config.gas.enableThermalMotion);
  }

  if (const auto* simulation = FindObjectSection(root, "simulation")) {
    config.simulation.numPrimaries =
        ReadOptionalSizeT(*simulation, "simulation", "num_primaries",
                          config.simulation.numPrimaries);
    config.simulation.initialElectronEnergyEv =
        ReadOptionalDouble(*simulation, "simulation", "initial_electron_energy_ev",
                           config.simulation.initialElectronEnergyEv);
    config.simulation.maxElectronEnergyEv =
        ReadOptionalDouble(*simulation, "simulation", "max_electron_energy_ev",
                           config.simulation.maxElectronEnergyEv);
    config.simulation.avalancheSizeLimit =
        ReadOptionalSizeT(*simulation, "simulation", "avalanche_size_limit",
                          config.simulation.avalancheSizeLimit);
    config.simulation.temperatureK =
        ReadOptionalDouble(*simulation, "simulation", "temperature_k",
                           config.simulation.temperatureK);
  }

  if (const auto* histogram = FindObjectSection(root, "histogram")) {
    config.histogram.energyBins =
        ReadOptionalInt(*histogram, "histogram", "energy_bins",
                        config.histogram.energyBins);
    config.histogram.multiplicityBinsMin =
        ReadOptionalInt(*histogram, "histogram", "multiplicity_bins_min",
                        config.histogram.multiplicityBinsMin);
  }

  if (config.baseline.distanceNm <= 0.) {
    throw std::runtime_error("baseline.distance_nm must be positive.");
  }
  if (config.baseline.pressureTorr <= 0.) {
    throw std::runtime_error("baseline.pressure_torr must be positive.");
  }
  if (config.baseline.fieldKvPerCm <= 0.) {
    throw std::runtime_error("baseline.field_kv_per_cm must be positive.");
  }
  if (config.scan.distanceNm.empty()) {
    throw std::runtime_error("scan.distance_nm must contain at least one value.");
  }
  if (config.scan.pressureTorr.empty()) {
    throw std::runtime_error("scan.pressure_torr must contain at least one value.");
  }
  if (config.scan.fieldKvPerCm.empty()) {
    throw std::runtime_error("scan.field_kv_per_cm must contain at least one value.");
  }
  if (config.simulation.numPrimaries == 0) {
    throw std::runtime_error("simulation.num_primaries must be at least 1.");
  }
  if (config.simulation.maxElectronEnergyEv <= 0.) {
    throw std::runtime_error("simulation.max_electron_energy_ev must be positive.");
  }
  if (config.simulation.initialElectronEnergyEv < 0.) {
    throw std::runtime_error(
        "simulation.initial_electron_energy_ev must be non-negative.");
  }
  if (config.simulation.temperatureK <= 0.) {
    throw std::runtime_error("simulation.temperature_k must be positive.");
  }
  if (config.histogram.energyBins <= 0) {
    throw std::runtime_error("histogram.energy_bins must be positive.");
  }
  if (config.histogram.multiplicityBinsMin <= 0) {
    throw std::runtime_error("histogram.multiplicity_bins_min must be positive.");
  }

  return config;
}

std::vector<ScanKind> ResolveScans(const std::string& scanSpec) {
  const std::string spec = ToLower(scanSpec);
  if (spec == "all") {
    return {ScanKind::Distance, ScanKind::Pressure, ScanKind::Field};
  }
  if (spec == "d" || spec == "distance") return {ScanKind::Distance};
  if (spec == "p" || spec == "pressure") return {ScanKind::Pressure};
  if (spec == "field" || spec == "e") return {ScanKind::Field};
  throw std::runtime_error(
      "Unknown scan selector '" + scanSpec + "'. Use all, d, p, or field.");
}

bool HasScan(const std::vector<ScanKind>& scanKinds, const ScanKind target) {
  return std::find(scanKinds.begin(), scanKinds.end(), target) != scanKinds.end();
}

std::string CanonicalScanLabel(const std::vector<ScanKind>& scanKinds) {
  if (scanKinds.size() == 3 && HasScan(scanKinds, ScanKind::Distance) &&
      HasScan(scanKinds, ScanKind::Pressure) && HasScan(scanKinds, ScanKind::Field)) {
    return "all";
  }
  if (scanKinds.size() != 1) {
    throw std::runtime_error("Unsupported scan combination when naming the run folder.");
  }
  switch (scanKinds.front()) {
    case ScanKind::Distance:
      return "distance";
    case ScanKind::Pressure:
      return "pressure";
    case ScanKind::Field:
      return "field";
  }
  return "scan";
}

std::string DistanceToken(const Config& config, const std::vector<ScanKind>& scanKinds) {
  if (HasScan(scanKinds, ScanKind::Distance)) {
    return JoinPathTokenValues(config.scan.distanceNm, "nm");
  }
  return FileSafeNumber(config.baseline.distanceNm) + "nm";
}

std::string PressureToken(const Config& config, const std::vector<ScanKind>& scanKinds) {
  if (HasScan(scanKinds, ScanKind::Pressure)) {
    return JoinPathTokenValues(config.scan.pressureTorr, "Torr");
  }
  return FileSafeNumber(config.baseline.pressureTorr) + "Torr";
}

std::string FieldToken(const Config& config, const std::vector<ScanKind>& scanKinds) {
  if (HasScan(scanKinds, ScanKind::Field)) {
    return JoinPathTokenValues(config.scan.fieldKvPerCm, "kVcm");
  }
  return FileSafeNumber(config.baseline.fieldKvPerCm) + "kVcm";
}

std::string BuildRunFolderName(const Config& config,
                               const std::vector<ScanKind>& scanKinds) {
  std::ostringstream stream;
  stream << "scan-" << CanonicalScanLabel(scanKinds)
         << "__d-" << DistanceToken(config, scanKinds)
         << "__p-" << PressureToken(config, scanKinds)
         << "__E-" << FieldToken(config, scanKinds)
         << "__n-" << config.simulation.numPrimaries;
  return stream.str();
}

ScanPoint MakeDistancePoint(const Config& config, const double distanceNm) {
  ScanPoint point;
  point.kind = ScanKind::Distance;
  point.scanKey = "distance";
  point.parameterKey = "distance_nm";
  point.parameterUnit = "nm";
  point.xAxisTitle = "Gap d [nm]";
  point.parameterValue = distanceNm;
  point.distanceNm = distanceNm;
  point.pressureTorr = config.baseline.pressureTorr;
  point.fieldKvPerCm = config.baseline.fieldKvPerCm;
  point.tag = "distance_d_" + FileSafeNumber(distanceNm) + "nm";
  point.label = "d = " + FormatNumber(distanceNm) + " nm";
  return point;
}

ScanPoint MakePressurePoint(const Config& config, const double pressureTorr) {
  ScanPoint point;
  point.kind = ScanKind::Pressure;
  point.scanKey = "pressure";
  point.parameterKey = "pressure_torr";
  point.parameterUnit = "Torr";
  point.xAxisTitle = "Pressure [Torr]";
  point.parameterValue = pressureTorr;
  point.distanceNm = config.baseline.distanceNm;
  point.pressureTorr = pressureTorr;
  point.fieldKvPerCm = config.baseline.fieldKvPerCm;
  point.tag = "pressure_p_" + FileSafeNumber(pressureTorr) + "Torr";
  point.label = "p = " + FormatNumber(pressureTorr) + " Torr";
  return point;
}

ScanPoint MakeFieldPoint(const Config& config, const double fieldKvPerCm) {
  ScanPoint point;
  point.kind = ScanKind::Field;
  point.scanKey = "field";
  point.parameterKey = "field_kv_per_cm";
  point.parameterUnit = "kV/cm";
  point.xAxisTitle = "Field [kV/cm]";
  point.parameterValue = fieldKvPerCm;
  point.distanceNm = config.baseline.distanceNm;
  point.pressureTorr = config.baseline.pressureTorr;
  point.fieldKvPerCm = fieldKvPerCm;
  point.tag = "field_E_" + FileSafeNumber(fieldKvPerCm) + "kVcm";
  point.label = "E = " + FormatNumber(fieldKvPerCm) + " kV/cm";
  return point;
}

std::vector<ScanPoint> BuildScanPoints(const Config& config, const ScanKind kind) {
  std::vector<ScanPoint> points;
  switch (kind) {
    case ScanKind::Distance:
      points.reserve(config.scan.distanceNm.size());
      for (const double value : config.scan.distanceNm) {
        points.push_back(MakeDistancePoint(config, value));
      }
      break;
    case ScanKind::Pressure:
      points.reserve(config.scan.pressureTorr.size());
      for (const double value : config.scan.pressureTorr) {
        points.push_back(MakePressurePoint(config, value));
      }
      break;
    case ScanKind::Field:
      points.reserve(config.scan.fieldKvPerCm.size());
      for (const double value : config.scan.fieldKvPerCm) {
        points.push_back(MakeFieldPoint(config, value));
      }
      break;
  }
  return points;
}

PointSummary SummarisePoint(const ScanPoint& point, const Config& config,
                           const std::vector<double>& energies,
                           const std::vector<double>& multiplicities,
                           const std::size_t primariesWithArrival,
                           const std::size_t sizeLimitHitPrimaries,
                           const std::size_t transportFailures) {
  PointSummary summary;
  summary.scanKey = point.scanKey;
  summary.parameterKey = point.parameterKey;
  summary.parameterUnit = point.parameterUnit;
  summary.parameterValue = point.parameterValue;
  summary.numPrimaries = config.simulation.numPrimaries;
  summary.totalArrivals = energies.size();
  summary.primariesWithArrival = primariesWithArrival;
  summary.sizeLimitHitPrimaries = sizeLimitHitPrimaries;
  summary.transportFailures = transportFailures;
  summary.arrivalFraction =
      summary.numPrimaries == 0
          ? 0.
          : static_cast<double>(summary.primariesWithArrival) /
                static_cast<double>(summary.numPrimaries);
  summary.meanEnergyEv = Mean(energies);
  summary.rmsEnergyEv = Rms(energies, summary.meanEnergyEv);
  summary.semEnergyEv = Sem(summary.rmsEnergyEv, energies.size());
  summary.meanMultiplicity = Mean(multiplicities);
  summary.rmsMultiplicity = Rms(multiplicities, summary.meanMultiplicity);
  summary.semMultiplicity = Sem(summary.rmsMultiplicity, multiplicities.size());
  summary.xAxisTitle = point.xAxisTitle;
  return summary;
}

TH1D MakeEnergyHistogram(const std::string& name, const Config& config,
                         const std::vector<double>& energies) {
  TH1D histogram(name.c_str(), "", config.histogram.energyBins, 0.,
                 config.simulation.maxElectronEnergyEv);
  histogram.SetDirectory(nullptr);
  histogram.SetTitle("Arrival electron energy;Final electron energy [eV];Count");
  for (const double energy : energies) histogram.Fill(energy);
  return histogram;
}

TH1D MakeMultiplicityHistogram(const std::string& name, const Config& config,
                               const std::vector<double>& multiplicities) {
  double maxMultiplicity = 0.;
  for (const double value : multiplicities) {
    if (value > maxMultiplicity) maxMultiplicity = value;
  }
  const int bins = std::max(config.histogram.multiplicityBinsMin,
                            static_cast<int>(std::ceil(maxMultiplicity)) + 1);
  TH1D histogram(name.c_str(), "", bins, -0.5, static_cast<double>(bins) - 0.5);
  histogram.SetDirectory(nullptr);
  histogram.SetTitle(
      "Arrival multiplicity per primary;Arriving electrons per primary;Count");
  for (const double multiplicity : multiplicities) histogram.Fill(multiplicity);
  return histogram;
}

void SaveHistogramPair(const TH1D& energyHistogram, const TH1D& multiplicityHistogram,
                       const ScanPoint& point, const fs::path& outputPath) {
  TCanvas canvas(("c_" + point.tag).c_str(), point.label.c_str(), 1200, 500);
  canvas.Divide(2, 1);

  canvas.cd(1);
  TH1D energyCopy = energyHistogram;
  energyCopy.SetTitle((point.label + ";Final electron energy [eV];Count").c_str());
  energyCopy.SetLineWidth(2);
  energyCopy.Draw("HIST");

  canvas.cd(2);
  TH1D multiplicityCopy = multiplicityHistogram;
  multiplicityCopy.SetTitle(
      (point.label + ";Arriving electrons per primary;Count").c_str());
  multiplicityCopy.SetLineWidth(2);
  multiplicityCopy.Draw("HIST");

  canvas.SaveAs(outputPath.string().c_str());
}

void WritePointHistograms(TDirectory* scanDir, const ScanPoint& point,
                          const TH1D& energyHistogram,
                          const TH1D& multiplicityHistogram) {
  TDirectory* pointDir = scanDir->mkdir(point.tag.c_str());
  if (!pointDir) {
    throw std::runtime_error("Failed to create ROOT directory for " + point.tag);
  }
  pointDir->cd();

  TH1D energyCopy = energyHistogram;
  energyCopy.SetName("arrival_energy");
  energyCopy.Write();

  TH1D multiplicityCopy = multiplicityHistogram;
  multiplicityCopy.SetName("arrival_multiplicity");
  multiplicityCopy.Write();
}

PointData RunScanPoint(const Config& config, const ScanPoint& point,
                       TDirectory* scanDir, const fs::path& histogramDir) {
  const double distanceCm = NmToCm(point.distanceNm);
  const double fieldVPerCm = KvPerCmToVPerCm(point.fieldKvPerCm);
  const double deltaVoltage = fieldVPerCm * distanceCm;
  const double sourceX = 0.01 * distanceCm;
  const double hitTolerance = std::max(1.e-8, 0.02 * distanceCm);
  const double lateralHalfWidth = std::max(100. * distanceCm, 1.e-3);

  MediumMagboltz gas("ar", 93., "co2", 7.);
  gas.SetTemperature(config.simulation.temperatureK);
  gas.SetPressure(point.pressureTorr);
  gas.SetMaxElectronEnergy(config.simulation.maxElectronEnergyEv);
  gas.EnableThermalMotion(config.gas.enableThermalMotion);

  if (!gas.Initialise(false)) {
    throw std::runtime_error("Failed to initialise MediumMagboltz for scan point " +
                             point.label);
  }
  if (config.gas.enablePenning) {
    if (!gas.EnablePenningTransfer()) {
      std::cerr << "Warning: Penning transfer could not be enabled for "
                << point.label << ". Continuing without it.\n";
    }
  } else {
    gas.DisablePenningTransfer();
  }

  ComponentAnalyticField component;
  component.SetMedium(&gas);
  component.AddPlaneX(0., 0., "source");
  component.AddPlaneX(distanceCm, deltaVoltage, "arrival");

  Sensor sensor(&component);
  if (!sensor.SetArea(-0.1 * distanceCm, -lateralHalfWidth, -lateralHalfWidth,
                      1.1 * distanceCm, lateralHalfWidth, lateralHalfWidth)) {
    throw std::runtime_error("Failed to set the sensor area for " + point.label);
  }

  AvalancheMicroscopic avalanche(&sensor);
  avalanche.EnableSignalCalculation(false);
  if (config.simulation.avalancheSizeLimit > 0) {
    avalanche.EnableAvalancheSizeLimit(config.simulation.avalancheSizeLimit);
  }

  std::vector<double> energies;
  std::vector<double> multiplicities;
  energies.reserve(config.simulation.numPrimaries);
  multiplicities.reserve(config.simulation.numPrimaries);

  std::size_t primariesWithArrival = 0;
  std::size_t sizeLimitHitPrimaries = 0;
  std::size_t transportFailures = 0;
  EndpointStatusBreakdown statusBreakdown;

  const std::size_t progressStep =
      std::max<std::size_t>(1, config.simulation.numPrimaries / 5);

  for (std::size_t primary = 0; primary < config.simulation.numPrimaries; ++primary) {
    const bool ok = avalanche.AvalancheElectron(
        sourceX, 0., 0., 0., config.simulation.initialElectronEnergyEv);
    if (!ok) ++transportFailures;

    int nElectrons = 0;
    int nIons = 0;
    avalanche.GetAvalancheSize(nElectrons, nIons);
    if (config.simulation.avalancheSizeLimit > 0 &&
        nElectrons >= static_cast<int>(config.simulation.avalancheSizeLimit)) {
      ++sizeLimitHitPrimaries;
    }

    std::size_t arrivalsThisPrimary = 0;
    const auto endpoints = avalanche.GetNumberOfElectronEndpoints();
    for (std::size_t i = 0; i < endpoints; ++i) {
      double x0 = 0., y0 = 0., z0 = 0., t0 = 0., e0 = 0.;
      double x1 = 0., y1 = 0., z1 = 0., t1 = 0., e1 = 0.;
      int status = 0;
      avalanche.GetElectronEndpoint(i, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1,
                                    status);
      const bool atArrivalBoundary = std::abs(x1 - distanceCm) <= hitTolerance;
      const bool atSourceBoundary = std::abs(x1) <= hitTolerance;
      const bool hitArrivalPlane =
          atArrivalBoundary &&
          (status == Garfield::StatusHitPlane ||
           status == Garfield::StatusLeftDriftMedium);
      if (status == Garfield::StatusHitPlane) {
        if (atArrivalBoundary) {
          ++statusBreakdown.hitArrivalPlane;
        } else if (atSourceBoundary) {
          ++statusBreakdown.hitSourcePlane;
        } else {
          ++statusBreakdown.hitOtherPlane;
        }
      } else if (status == Garfield::StatusAttached) {
        ++statusBreakdown.attached;
      } else if (status == Garfield::StatusBelowTransportCut) {
        ++statusBreakdown.belowTransportCut;
      } else if (status == Garfield::StatusLeftDriftMedium) {
        if (atArrivalBoundary) {
          ++statusBreakdown.leftDriftMediumArrivalBoundary;
        } else if (atSourceBoundary) {
          ++statusBreakdown.leftDriftMediumSourceBoundary;
        } else {
          ++statusBreakdown.leftDriftMediumOther;
        }
      } else if (status == Garfield::StatusLeftDriftArea) {
        ++statusBreakdown.leftDriftArea;
      } else if (status == Garfield::StatusOutsideTimeWindow) {
        ++statusBreakdown.outsideTimeWindow;
      } else if (status == Garfield::StatusCalculationAbandoned) {
        ++statusBreakdown.calculationAbandoned;
      } else {
        ++statusBreakdown.other;
      }
      if (!hitArrivalPlane) continue;
      energies.push_back(e1);
      ++arrivalsThisPrimary;
    }

    multiplicities.push_back(static_cast<double>(arrivalsThisPrimary));
    if (arrivalsThisPrimary > 0) ++primariesWithArrival;

    if ((primary + 1) % progressStep == 0 ||
        primary + 1 == config.simulation.numPrimaries) {
      std::cout << "  " << point.label << ": " << (primary + 1) << "/"
                << config.simulation.numPrimaries << " primaries processed\n";
    }
  }

  PointData pointData;
  pointData.arrivalEnergiesEv = std::move(energies);
  pointData.arrivalMultiplicities = std::move(multiplicities);
  pointData.summary =
      SummarisePoint(point, config, pointData.arrivalEnergiesEv,
                     pointData.arrivalMultiplicities, primariesWithArrival,
                     sizeLimitHitPrimaries, transportFailures);

  TH1D energyHistogram =
      MakeEnergyHistogram("arrival_energy", config, pointData.arrivalEnergiesEv);
  TH1D multiplicityHistogram = MakeMultiplicityHistogram(
      "arrival_multiplicity", config, pointData.arrivalMultiplicities);

  WritePointHistograms(scanDir, point, energyHistogram, multiplicityHistogram);

  EnsureDirectory(histogramDir);
  SaveHistogramPair(energyHistogram, multiplicityHistogram, point,
                    histogramDir / (point.tag + ".png"));

  std::cout << "    arrivals: " << pointData.summary.totalArrivals
            << ", mean energy: " << FormatNumber(pointData.summary.meanEnergyEv)
            << " eV, mean multiplicity: "
            << FormatNumber(pointData.summary.meanMultiplicity) << "\n";
  if (pointData.summary.sizeLimitHitPrimaries > 0) {
    std::cout << "    size limit hit in "
              << pointData.summary.sizeLimitHitPrimaries << " / "
              << pointData.summary.numPrimaries << " primaries\n";
  }
  if (pointData.summary.transportFailures > 0) {
    std::cout << "    transport returned failure in "
              << pointData.summary.transportFailures << " primaries\n";
  }
  std::cout << "    endpoint statuses:"
            << " hit-arrival=" << statusBreakdown.hitArrivalPlane
            << ", exit-arrival=" << statusBreakdown.leftDriftMediumArrivalBoundary
            << ", hit-source=" << statusBreakdown.hitSourcePlane
            << ", exit-source=" << statusBreakdown.leftDriftMediumSourceBoundary
            << ", attached=" << statusBreakdown.attached
            << ", below-cut=" << statusBreakdown.belowTransportCut
            << ", left-area=" << statusBreakdown.leftDriftArea
            << ", abandoned=" << statusBreakdown.calculationAbandoned
            << ", other=" << statusBreakdown.other + statusBreakdown.hitOtherPlane +
                                   statusBreakdown.leftDriftMediumOther
            << "\n";

  return pointData;
}

void WriteSummaryGraphs(const std::vector<PointSummary>& summaries, TDirectory* scanDir,
                        const fs::path& outputPath) {
  if (summaries.empty()) return;

  std::vector<double> x(summaries.size());
  std::vector<double> xerr(summaries.size(), 0.);
  std::vector<double> energy(summaries.size());
  std::vector<double> energyErr(summaries.size());
  std::vector<double> multiplicity(summaries.size());
  std::vector<double> multiplicityErr(summaries.size());

  for (std::size_t i = 0; i < summaries.size(); ++i) {
    x[i] = summaries[i].parameterValue;
    energy[i] = summaries[i].meanEnergyEv;
    energyErr[i] = summaries[i].semEnergyEv;
    multiplicity[i] = summaries[i].meanMultiplicity;
    multiplicityErr[i] = summaries[i].semMultiplicity;
  }

  const std::string scanKey = summaries.front().scanKey;
  const std::string xAxisTitle = summaries.front().xAxisTitle;

  TGraphErrors energyGraph(static_cast<int>(summaries.size()), x.data(),
                           energy.data(), xerr.data(), energyErr.data());
  energyGraph.SetName(("mean_arrival_energy_vs_" + scanKey).c_str());
  energyGraph.SetTitle(
      ("Mean arrival energy;" + xAxisTitle + ";Mean arrival energy [eV]").c_str());
  energyGraph.SetMarkerStyle(20);
  energyGraph.SetLineWidth(2);

  TGraphErrors multiplicityGraph(static_cast<int>(summaries.size()), x.data(),
                                 multiplicity.data(), xerr.data(),
                                 multiplicityErr.data());
  multiplicityGraph.SetName(("mean_arrival_multiplicity_vs_" + scanKey).c_str());
  multiplicityGraph.SetTitle(("Mean arrival multiplicity;" + xAxisTitle +
                              ";Mean arrival multiplicity")
                                 .c_str());
  multiplicityGraph.SetMarkerStyle(21);
  multiplicityGraph.SetLineWidth(2);

  TDirectory* summaryDir = scanDir->mkdir("summary");
  if (!summaryDir) {
    throw std::runtime_error("Failed to create ROOT summary directory for " +
                             scanKey);
  }
  summaryDir->cd();
  energyGraph.Write();
  multiplicityGraph.Write();

  TCanvas canvas(("c_summary_" + scanKey).c_str(), scanKey.c_str(), 1200, 500);
  canvas.Divide(2, 1);
  canvas.cd(1);
  energyGraph.Draw("APL");
  canvas.cd(2);
  multiplicityGraph.Draw("APL");
  canvas.SaveAs(outputPath.string().c_str());
}

void WriteSummaryCsv(const fs::path& csvPath, const std::vector<PointSummary>& summaries) {
  std::ofstream stream(csvPath);
  if (!stream) {
    throw std::runtime_error("Failed to open CSV file for writing: " +
                             csvPath.string());
  }

  stream << "scan,parameter,parameter_unit,parameter_value,num_primaries,"
            "total_arrivals,primaries_with_arrival,arrival_fraction,"
            "mean_arrival_energy_ev,rms_arrival_energy_ev,sem_arrival_energy_ev,"
            "mean_multiplicity,rms_multiplicity,sem_multiplicity,"
            "size_limit_hit_primaries,transport_failures\n";

  stream << std::fixed << std::setprecision(6);
  for (const auto& summary : summaries) {
    stream << summary.scanKey << ','
           << summary.parameterKey << ','
           << summary.parameterUnit << ','
           << summary.parameterValue << ','
           << summary.numPrimaries << ','
           << summary.totalArrivals << ','
           << summary.primariesWithArrival << ','
           << summary.arrivalFraction << ','
           << summary.meanEnergyEv << ','
           << summary.rmsEnergyEv << ','
           << summary.semEnergyEv << ','
           << summary.meanMultiplicity << ','
           << summary.rmsMultiplicity << ','
           << summary.semMultiplicity << ','
           << summary.sizeLimitHitPrimaries << ','
           << summary.transportFailures << '\n';
  }
}

std::string ScanTitle(const ScanKind kind) {
  switch (kind) {
    case ScanKind::Distance:
      return "distance";
    case ScanKind::Pressure:
      return "pressure";
    case ScanKind::Field:
      return "field";
  }
  return "scan";
}

json ConfigToJson(const Config& config) {
  return json{
      {"baseline",
       {{"distance_nm", config.baseline.distanceNm},
        {"pressure_torr", config.baseline.pressureTorr},
        {"field_kv_per_cm", config.baseline.fieldKvPerCm}}},
      {"scan",
       {{"distance_nm", config.scan.distanceNm},
        {"pressure_torr", config.scan.pressureTorr},
        {"field_kv_per_cm", config.scan.fieldKvPerCm}}},
      {"gas",
       {{"enable_penning", config.gas.enablePenning},
        {"enable_thermal_motion", config.gas.enableThermalMotion}}},
      {"simulation",
       {{"num_primaries", config.simulation.numPrimaries},
        {"initial_electron_energy_ev",
         config.simulation.initialElectronEnergyEv},
        {"max_electron_energy_ev", config.simulation.maxElectronEnergyEv},
        {"avalanche_size_limit", config.simulation.avalancheSizeLimit},
        {"temperature_k", config.simulation.temperatureK}}},
      {"histogram",
       {{"energy_bins", config.histogram.energyBins},
        {"multiplicity_bins_min", config.histogram.multiplicityBinsMin}}}};
}

void WriteJsonFile(const fs::path& path, const json& payload) {
  std::ofstream stream(path);
  if (!stream) {
    throw std::runtime_error("Failed to open JSON file for writing: " +
                             path.string());
  }
  stream << std::setw(2) << payload << '\n';
}

void WriteConfigUsed(const fs::path& outputDir, const Config& config) {
  WriteJsonFile(outputDir / "run_config_used.json", ConfigToJson(config));
}

void WriteRunManifest(const fs::path& metadataPath, const fs::path& configPath,
                      const fs::path& outputBaseDir, const fs::path& runOutputDir,
                      const std::vector<ScanKind>& scanKinds,
                      const std::string& commandLine, const Config& config) {
  const json manifest = {
      {"scan", CanonicalScanLabel(scanKinds)},
      {"run_folder", runOutputDir.filename().string()},
      {"output_base_dir", outputBaseDir.string()},
      {"run_output_dir", runOutputDir.string()},
      {"config_path", configPath.string()},
      {"num_primaries", config.simulation.numPrimaries},
      {"command_line", commandLine},
      {"baseline",
       {{"distance_nm", config.baseline.distanceNm},
        {"pressure_torr", config.baseline.pressureTorr},
        {"field_kv_per_cm", config.baseline.fieldKvPerCm}}},
      {"scan_grid",
       {{"distance_nm", config.scan.distanceNm},
        {"pressure_torr", config.scan.pressureTorr},
        {"field_kv_per_cm", config.scan.fieldKvPerCm}}},
      {"artifacts",
       {{"root_file", "parallel_plate_scan.root"},
        {"summary_csv", "scan_summary.csv"},
        {"config_used_json", "run_config_used.json"},
        {"manifest_json", "run_manifest.json"},
        {"summary_dir", "summary"},
        {"histograms_dir", "histograms"}}}};
  WriteJsonFile(metadataPath, manifest);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(1110);
    TH1::AddDirectory(false);
    TH1::StatOverflows(true);
    gRandom->SetSeed(0);

    const CliOptions options = ParseCommandLine(argc, argv);
    const Config config = LoadConfig(options.configPath);
    const std::vector<ScanKind> selectedScans = ResolveScans(options.scanSpec);
    const std::string commandLine = BuildCommandLine(argc, argv);
    const fs::path runOutputDir =
        options.outDir / BuildRunFolderName(config, selectedScans);

    EnsureDirectory(options.outDir);
    EnsureDirectory(runOutputDir);
    EnsureDirectory(runOutputDir / "summary");
    EnsureDirectory(runOutputDir / "histograms");
    WriteConfigUsed(runOutputDir, config);
    WriteRunManifest(runOutputDir / "run_manifest.json", options.configPath,
                     options.outDir, runOutputDir, selectedScans, commandLine,
                     config);

    TFile outputFile((runOutputDir / "parallel_plate_scan.root").string().c_str(),
                     "RECREATE");
    if (outputFile.IsZombie()) {
      throw std::runtime_error("Failed to create ROOT output file in " +
                               runOutputDir.string());
    }

    std::vector<PointSummary> allSummaries;

    std::cout << "Running Garfield++ nanoscale parallel-plate scan\n";
    std::cout << "  config: " << options.configPath << "\n";
    std::cout << "  output base: " << options.outDir << "\n";
    std::cout << "  run folder: " << runOutputDir << "\n";

    for (const ScanKind scanKind : selectedScans) {
      const auto points = BuildScanPoints(config, scanKind);
      const std::string scanKey = ScanTitle(scanKind);
      std::cout << "Starting " << scanKey << " scan with " << points.size()
                << " point(s)\n";

      TDirectory* scanDir = outputFile.mkdir((scanKey + "_scan").c_str());
      if (!scanDir) {
        throw std::runtime_error("Failed to create ROOT directory for scan " +
                                 scanKey);
      }

      std::vector<PointSummary> scanSummaries;
      scanSummaries.reserve(points.size());
      for (const auto& point : points) {
        PointData pointData = RunScanPoint(
            config, point, scanDir, runOutputDir / "histograms" / scanKey);
        scanSummaries.push_back(pointData.summary);
        allSummaries.push_back(pointData.summary);
      }

      WriteSummaryGraphs(scanSummaries, scanDir,
                         runOutputDir / "summary" / (scanKey + "_summary.png"));
    }

    outputFile.Write();
    outputFile.Close();

    WriteSummaryCsv(runOutputDir / "scan_summary.csv", allSummaries);

    std::cout << "Finished. Results written to " << runOutputDir << "\n";
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << "Error: " << exception.what() << "\n";
    return 1;
  }
}
