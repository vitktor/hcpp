#pragma once

namespace hc
{

enum class PathStatus
{
  Success,
  Diverged,
  Failed,
};

inline const char* to_string(PathStatus status)
{
  switch (status)
  {
  case PathStatus::Success:
    return "success";
  case PathStatus::Diverged:
    return "diverged";
  case PathStatus::Failed:
    return "failed";
  }
  return "unknown";
}

} // namespace hc
