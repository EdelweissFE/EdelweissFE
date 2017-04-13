#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 09:27:14 2017

@author: matthias
"""

from fe.outputmanagers.nodemonitor import OutputManager as nodeMonitor
from fe.outputmanagers.nodesetmonitor import OutputManager as nodeSetMonitor
from fe.outputmanagers.ensight import OutputManager as ensight
from fe.outputmanagers.meshplot import OutputManager as meshPlot
from fe.outputmanagers.timemonitor import OutputManager as timeMonitor

outputManagersLibrary = {
                    'nodemonitor' : nodeMonitor,
                    'nodesetmonitor': nodeSetMonitor,
                    'ensight': ensight,
                    'meshplot': meshPlot,
                    'timemonitor': timeMonitor
                    }
