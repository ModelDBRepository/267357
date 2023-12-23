: Crude threshold point process to allow for multiple 'netcon' like objects in a single compartment


NEURON {
  POINT_PROCESS SecNet
  RANGE thresh
}
PARAMETER { thresh = -60 } : threshold is -20 mV
ASSIGNED { v }
INITIAL {
  net_send(0, 1) : to execute the WATCH statement
}
NET_RECEIVE(w) {
  if (flag == 1) {
    WATCH (v > thresh) 2 : if v crosses thresh in a positive-going direction, send a self event with weight 2
  } else if (flag == 2) {
    net_event(t) : send an event to all NetCons that have this point process as their source
  }
}
