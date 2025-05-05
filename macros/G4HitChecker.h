#ifndef G4HITCHECKER_H
#define G4HITCHECKER_H

// Standard sPHENIX / Fun4All includes
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// G4Hit types
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

// Additional includes if we want run/event:
#include <ffaobjects/EventHeader.h>
#include <TSystem.h>

// C++ includes
#include <string>
#include <iostream>

/**
 * G4HitChecker
 *
 * A SubsysReco that checks for G4Hits in a specified detector node,
 * e.g. "CEMC" => node name "G4HIT_CEMC".
 * 
 * This class can also enforce an optional "maxAllowedZeroHitEvents" threshold:
 *  - If the number of events with zero hits in a row exceeds that threshold,
 *    the code can choose to gracefully exit or throw an error.
 */
class G4HitChecker : public SubsysReco
{
public:
  /**
   * Constructor
   * @param detName: short name of the detector, e.g. "CEMC".
   * @param maxZeroAllowed: if >0, the number of consecutive zero-hit events
   *                        allowed before we optionally bail out.
   */
  G4HitChecker(const std::string &detName = "CEMC",
               int maxZeroAllowed = 0)
    : SubsysReco("G4HitChecker")
    , m_detName(detName)
    , m_maxConsecutiveZeros(maxZeroAllowed)
    , m_currentZeroCount(0)
  {
    // Example: by default, set some moderate verbosity
    Verbosity(1);
  }

  /// SubsysReco method: Called once per event
  int process_event(PHCompositeNode *topNode) override
  {
    // Optionally retrieve run/event from EventHeader
    EventHeader* evtheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
    int thisEvent = (evtheader ? evtheader->get_EvtSequence() : -9999);
    int thisRun   = (evtheader ? evtheader->get_RunNumber()   : -9999);

    // Build the G4HIT node name, e.g. "G4HIT_CEMC"
    std::string nodeName = "G4HIT_" + m_detName;

    // Retrieve the hit container
    PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodeName.c_str());
    if (!hits)
    {
      if (Verbosity() > 0)
      {
        std::cout << "[G4HitChecker] WARNING: No PHG4HitContainer named '"
                  << nodeName << "' found on run=" << thisRun
                  << ", evt=" << thisEvent
                  << ". Possibly a geometry/input mismatch, or no hits exist.\n";
      }
      // This counts as zero hits
      handleZeroHits();
      return 0;
    }

    // How many hits in this container?
    unsigned int nhits = hits->size();

    // Basic debug print
    if (Verbosity() > 0)
    {
      std::cout << "[G4HitChecker] run=" << thisRun << ", evt=" << thisEvent
                << ": " << nodeName
                << " has " << nhits << " hits on this event." << std::endl;
    }

    // If advanced debugging is on, we might print out a few hits
    if (Verbosity() > 2 && nhits > 0)
    {
      std::cout << "   [G4HitChecker] Example listing of first 3 hits:\n";
      int count = 0;
      for (auto hitpair = hits->getHits().first;
                hitpair != hits->getHits().second; ++hitpair)
      {
        PHG4Hit *thishit = hitpair->second;
        std::cout << "     - hit id=" << thishit->get_hit_id()
                  << ", edep=" << thishit->get_edep()
                  << ", (x,y,z)=(" << thishit->get_avg_x() << ", "
                                   << thishit->get_avg_y() << ", "
                                   << thishit->get_avg_z() << ")"
                  << std::endl;
        count++;
        if(count >= 3) break;
      }
    }

    // If zero hits, handle logic
    if (nhits == 0)
    {
      handleZeroHits();
    }
    else
    {
      // reset the zero count
      m_currentZeroCount = 0;
    }

    return 0; // EVENT_OK
  }

private:
  /// Helper to handle 0-hit events, possibly counting them and exiting if threshold is exceeded.
  void handleZeroHits()
  {
    m_currentZeroCount++;
    if(m_maxConsecutiveZeros > 0 && m_currentZeroCount >= m_maxConsecutiveZeros)
    {
      std::cerr << "[G4HitChecker] ERROR: We have reached "
                << m_currentZeroCount << " consecutive events with zero hits "
                << "in G4HIT_" << m_detName << ". Exiting run." << std::endl;
      // Option 1: stop run gracefully:
      gSystem->Exit(1);

      // Option 2: throw an exception instead
      // throw std::runtime_error("Too many zero-hit events. Stopping.");
    }
  }

private:
  std::string m_detName;
  int m_maxConsecutiveZeros;
  int m_currentZeroCount;
};

#endif // G4HITCHECKER_H

