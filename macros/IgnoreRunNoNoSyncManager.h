// IgnoreRunNoNoSyncManager.h
#ifndef IGNORE_RUNNO_NOSYNC_MANAGER_H
#define IGNORE_RUNNO_NOSYNC_MANAGER_H

#include <fun4all/Fun4AllNoSyncDstInputManager.h>

// A custom input manager that explicitly ignores run number mismatch
class IgnoreRunNoNoSyncManager : public Fun4AllNoSyncDstInputManager
{
public:
  IgnoreRunNoNoSyncManager(const std::string &name = "NOSYNC_WITH_IGNORED_RUNNO",
                           const std::string &nodename = "TOP",
                           const std::string &topnodename = "TOP")
    : Fun4AllNoSyncDstInputManager(name, nodename, topnodename)
  {
  }

  // We override RunNumber() to ALWAYS return 0, so the Fun4AllSyncManager
  // sees "run=0" => no mismatch check is triggered.
  int RunNumber() const override
  {
    return 0; 
  }
};

#endif

