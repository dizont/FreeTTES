"""
1D stratified hot-water thermal energy storage (TES) tank model FreeTTES.

Original model developed by: Andreas Herwig 
Implemented in Python and further development by: 
Bogdan Narusavicius (bogdan.narusavicius@tu-dresden.de) @ TU Dresden GEWV

The core public entry point is :func:`main`.
"""

# // Code für Modellierung eines Großwärmespeichers

# // Import von benötigten Bibliotheken
import warnings
import os
#import sys
import json
from scipy import special, interpolate
from math import log, pi, tan, exp
import csv
# fuer excel datei lesen
import numpy as np


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# v2 split: config + I/O are isolated into dedicated modules
# ---------------------------------------------------------------------------

import FreeTTES_config as _cfg
import FreeTTES_io as _io

import logging
logger = logging.getLogger(__name__)

# Legacy-compatible globals expected by the legacy physics code
script_dir = _cfg.SCRIPT_DIR
abs_folder_temperaturProfile = _cfg.abs_folder_temperaturProfile
abs_outputs_pfad = _cfg.abs_outputs_pfad
sz_folder = _cfg.sz_folder
config_path = _cfg.config_path
speicher_param = _cfg.speicher_param

# Backwards compatible aliases for the moved I/O functions
__Modell_Ausgabe_Zeitschritt = _io.ausgabe_zeitschritt
__Modell_letzter_Zustand = _io.letzter_zustand
# // Speicherparameter definieren
Speicherzustand_ = {}
Fundamentzustand_ = {}
Kapazitaeten_ = {}
g = 9.81


# Sync legacy path globals from config module (v2 split)
def _sync_legacy_globals():
    global abs_folder_temperaturProfile, abs_outputs_pfad, sz_folder, config_path, speicher_param
    abs_folder_temperaturProfile = _cfg.abs_folder_temperaturProfile
    abs_outputs_pfad = _cfg.abs_outputs_pfad
    sz_folder = _cfg.sz_folder
    config_path = _cfg.config_path
    speicher_param = _cfg.speicher_param

vorgang = ""
flag = False

# // Routine, die die Speicherberechnungen durchführt und über ein Skript aufgerufen wird
def main(t, dt, m_VL,  m_RL, T_Zustrom, T_amb, eingabe_volumen=False, zustand_uebernehmen=False, zustand={}, T_Abstrom=0, T_DR=None, T_RL = 60):
    """
    Hauptroutine, in der die gesamte Speichersimulation für einen Zeitschritt stattfindet. Wird über ein separates Skript aufgerufen.

    Parameters:
        t (float): current simulated time (starts with 0) in hours
        dt (int): timestep in seconds
        m_VL (float): Massflow supply line in kg/s or m3/h depending on eingabe_volumen
        m_RL (float): Massflow return line (!=-m_VL)
        T_Zustrom (float): Temeperature of entering water in °C
        T_amb (float): Ambient temperature in °C
        eingabe_volumen (bool): Whether volumetric or massflow is given 
        zustand_uebernehmen (bool): Start simulation with a given state (measured)
        zustand (dict): the actual height:temperature measured values if zustand_uebernehmen
        T_Abstrom (float): for validation purpose only (or when Volumetric flow is given) in °C
        T_DR (float): Temperature of the steam layer above water, in case its dynamic in °C
        T_RL (float): Temperature which defines the minimal useful temperature in °C

    Returns:
        outputs: Dictionary will various output values
    """
    # // Vorlauf- und Rücklaufvolumenströme korrekt zuordnen

    # Normalerweise T_Abstrom ist eine gesuchte Größe, aber wenn man mit den Messwerten vergleicht,
    # man weisst vorher die T_Abstrom aus den Begleitdaten des Betreibers
    # die Eingaben des Betreibers sind nrmlws volumenstrom, und die das ganze Modell arbeitet mit Massestrom
    # also an dieser Stelle muss man die Größen berechnen. 
    if m_VL > 0:
        vorgang = "beladen"
        if eingabe_volumen:                                    
            m_VL = m_VL * __Modell_Stoffwerte("rho",T_Zustrom) 
            m_RL = m_RL * __Modell_Stoffwerte("rho",T_Abstrom) 
    elif m_VL < 0:
        vorgang = "entladen"
        if eingabe_volumen:
            m_VL = m_VL * __Modell_Stoffwerte("rho",T_Abstrom) 
            m_RL = m_RL * __Modell_Stoffwerte("rho",T_Zustrom)
    else:
        vorgang = "stillstand"
    logger.debug("vorgang=%s m_VL=%s m_RL=%s T_Zustrom=%s", vorgang, m_VL, m_RL, T_Zustrom)

    t = t
    dt = dt
    m_VL = m_VL
    m_RL = m_RL
    T_Zustrom = T_Zustrom
    T_amb = T_amb
    T_Abstrom = T_Abstrom

    # // initialize (config + output folders)
    _cfg.ensure_initialized(force_config_reload=(t == 0), T_DR=T_DR)
    # Config loading re-binds the underlying dict; resync module-level aliases.
    _sync_legacy_globals()
    if T_DR is not None:                                                       # falls DR-Temperatur vorgegeben wurde, wird diese in eine Variable geschrieben
        speicher_param["T_DR"] = T_DR

        

    # // Prüfen, ob Speichereintrittstemperatur innerhalb der Temperaturgrenzen liegt
    if T_Zustrom > 105 or T_Zustrom < 25:
        logger.error("T_Zustrom out of bounds: %s", T_Zustrom)
        raise ValueError(f"Speichereintrittstemperatur ({T_Zustrom}) liegt nicht in den Temperaturgrenzen zwischen 25 und 105 °C")

    # // Weitere Speicherparameter definieren
    # Diese Sachen hier finden keine Verwendung, aber können nutzvoll sein für tiefere Analyse
    Ausgabewerte = {}
    BilanzMasse = {}
    BilanzEnergie = 0
    Ausgabewerte["FuellStand"] = {}
    Ausgabewerte["p_unten"] = {}
    Ausgabewerte["m_Abw_global_in_kg"] = {}
    Ausgabewerte["H_Abw_global_in_kJ"] = {}
    Ausgabewerte["m_Abw_global_in_kg"] = {}
    Ausgabewerte["H_Abw_global_in_kJ"] = {}
    Ausgabewerte["check"] = {}

    macro_energie = 0

    Masse_Global_Modellgrenzen = 0
    Energie_Global_Modellgrenzen = 0
    E_Verlust_Mantel_alle_dt_sub = 0
    Q_oben = 0
    H_Abstrom = 0
    Fundamentzustand = {}
    Kapazitaeten = {}
    alle_Temperaturprofile  = {}

    # // Initialisierung des Speichers durchführen
    # \\ bei t = 0 mit vorgegebenen Speicherparametern
    if t == 0:
        initial_Zustand = __Modell_Initialisierung(
                                                speicher_param["H_UEB_start"],
                                                speicher_param["Beladefaktor_start"],
                                                0.00,
                                                alle_Temperaturprofile={}, spline_uebernehmen=zustand_uebernehmen, zustand=zustand)
        alle_Temperaturprofile = initial_Zustand[0]
        Speicherzustand = initial_Zustand[1]
        Fundamentzustand = initial_Zustand[2]
        Kapazitaeten = initial_Zustand[3]
    # \\ bei t > 0 aus letztem Zeitschritt
    # Frage: bei erneuter Abfrage von TRNSYS wird letzter gespeicherter Zustand aus der Berechnung zuvor genommen?
    else:
        Speicherzustand, Fundamentzustand, Kapazitaeten = __Modell_letzter_Zustand()
        h_WS = Speicherzustand[max(sorted(Speicherzustand))][1]/2 + max(sorted(Speicherzustand))


    # // Masse und Energie im Startzustand bestimmen
    # -------------------------------------------
    dh_mal_rho = [v[1] * __Modell_Stoffwerte("rho",v[0])                                # Produkt aus Höhe der Zelle und Dichte der Zelle = Masse bezogen auf die Querschnittsfläche
                 for k,v in sorted(Speicherzustand.items())]
    alle_H = [v[1] * __Modell_Stoffwerte("rho", v[0]) * __Modell_Stoffwerte("h",v[0])   # Enthalpie bezogen auf die Querschnittsfläche
             for k,v in sorted(Speicherzustand.items())]
    alle_C_fundament = [v[0] * v[1]                                                     # Kapazitäten bezogen auf die Querschnittsfläche im Fundament
                       * __Modell_Stoffwerte("rho_Fundament") 
                       * __Modell_Stoffwerte("cp_Fundament")  
                       for k,v in sorted(Fundamentzustand.items())]
    alle_C_speicher = [v[0] * v[2] for k,v in sorted(Kapazitaeten.items())]             # Kapazitäten bezogen auf die Querschnittsfläche im Speicher # Frage: Was steht in v[2]?

    # Masse und Energie im Startzustand bestimmen
    Masse_Global_Speicherzustand = speicher_param["A_Quer"] * sum(dh_mal_rho)
    Energie_Global_Speicherzustand = speicher_param["A_Quer"] * sum(alle_H)
    Energie_Global_Speicherzustand += speicher_param["A_Quer"] * sum(alle_C_fundament)
    Energie_Global_Speicherzustand += sum(alle_C_speicher)

    # // Schleife ueber alle Sub-Zeitschritte
    # ------------------------------------
    dt_Modell = 60                                          # angestrebte Subzeitschrittweite
    n_sub = round(dt / dt_Modell)                           # Anzahl an Subzeitschritten - berechnet mit Zeitschrittweite des Modells und übergebenem Zeitschritt von außen
    dt_sub = dt / n_sub                                     # finale Subzeitschrittweite
    for j in range(1, n_sub+1):                             # hier startet ein Subzeitschritt
        counterInv = 0
        aktuellSekunden = t * 3600
        ausgabezeit = (aktuellSekunden + dt_sub * j) / 3600 # in perl es ist string
        #print(j, ausgabezeit)

        # // Fall: Entladung durchführen - Zustrom

        if m_RL > 0: # Massestrom im Rücklauf ist positiv
            Begleitdaten = {}
            Begleitdaten[aktuellSekunden] = {}
            Begleitdaten[aktuellSekunden]["T_RL"] = T_Zustrom                       # T_RL wird die Zustromtemp. zugewiesen
            Begleitdaten[aktuellSekunden]["m_Punkt"] = m_RL
            Begleitdaten[aktuellSekunden]["m_Punkt_abstrom"] = m_VL

            # \\ Zustrommodell durchführen
            Speicherzustand = __Modell_Zustrom("unten", aktuellSekunden, ausgabezeit,
                                               Begleitdaten, dt_sub, Ausgabewerte,
                                               Speicherzustand)                     # neuen Speicherzustand nach Hinzufügen von Zellen unten berechnen
            Speicherzustand = __Modell_Aufraumen(Speicherzustand)                   # manchmal noetig weil dh wird zu 0
            # \\ Inversionen auflösen
            counterInv = 0
            Vp_zu = m_RL / __Modell_Stoffwerte("rho", T_Zustrom)
            start_inversion_status = __Modell_Inversionspruefung(Speicherzustand)   # überprüfen, ob Inversionen vorhanden sind
            inversion_status = start_inversion_status
            while inversion_status != "keine":                                      # so lange Inversionen vorhanden sind, werden diese in der Schleife aufgelöst
                counterInv += 1
                Speicherzustand = __Modell_Inversion(start_inversion_status,
                                                     "unten", 
                                                     Vp_zu, 
                                                     dt_sub, 
                                                     Speicherzustand)               #  Inversionen werden aufgelöst
                Speicherzustand = __Modell_Aufraumen(Speicherzustand)               # Aufräumen des Speichers
                if counterInv == 1:                                                 # beim ersten Durchlauf der Schleife wird das Impulsmodell durchgeführt
                    Speicherzustand = __Modell_Impuls(start_inversion_status,
                                                      dt_sub,
                                                      Speicherzustand)
                    # B: ACHTUNG es wird innerhalb impuls aufgerauemt, wieso nochmal?
                    Speicherzustand = __Modell_Aufraumen(Speicherzustand)           # Aufräumen des Speichers
                inversion_status = __Modell_Inversionspruefung(Speicherzustand)     # erneute Prüfung auf Inversionen
            #ende while

            # \\ Horizontalmischung durchführen
            Speicherzustand = __Modell_Horizontalmischung(Speicherzustand)

            for v in Speicherzustand.values():                                      # Impuls auf null setzen falls es noch nicht passiert ist
                v[2] = 0
                v[3] = 0

            # \\ Masse und Energie im Speicher neu berechnen
            masse_zu = m_RL * dt_sub
            Masse_Global_Modellgrenzen += masse_zu
            Energie_Global_Modellgrenzen += masse_zu\
                                            * __Modell_Stoffwerte("h", T_Zustrom)
        # ende if m_RL > 0

        # // Fall: Beladung durchführen - Abstrom

        elif m_RL < 0: # Massestrom im Rücklauf ist negativ
            Begleitdaten = {}
            Begleitdaten[aktuellSekunden] = {}
            Begleitdaten[aktuellSekunden]["m_Punkt"] = -m_RL                # Betrag
            Begleitdaten[aktuellSekunden]["T_Abstrom"] = T_Abstrom

            # \\ Abstrommodell durchführen
            Speicherzustand = __Modell_Abstrom("unten", aktuellSekunden,
                                               ausgabezeit, Ausgabewerte,
                                               Begleitdaten, dt_sub,
                                               Speicherzustand)             # neuen Speicherzustand nach Verkleinern von Zellen unten berechnen

            # \\ Masse und Energie im Speicher neu berechnen
            masse_ab = m_RL * dt_sub
            # B: ACHTUNG: noch nicht klar ob die Ausgabewerte so aussehen werden!
            theta_ab = Ausgabewerte["T_RL"][ausgabezeit]
            #theta_ab = T_Abstrom
            h_ab = __Modell_Stoffwerte("h", theta_ab)
            Masse_Global_Modellgrenzen += masse_ab    
            Energie_Global_Modellgrenzen += masse_ab * h_ab
            H_Abstrom += -masse_ab * h_ab
            # zugefuehrte energie
            micro_energie = __Modell_Stoffwerte("cp", (T_Zustrom+theta_ab)/2 ) * m_VL * (T_Zustrom - theta_ab) * dt_sub
            macro_energie += micro_energie
        # ende if m_RL<0

        # // Fall: Beladung durchführen - Zustrom

        if m_VL > 0: # Massestrom im Vorlauf ist positiv
            Begleitdaten = {}
            Begleitdaten[aktuellSekunden] = {}
            Begleitdaten[aktuellSekunden]["T_VL"] = T_Zustrom
            Begleitdaten[aktuellSekunden]["m_Punkt"] = m_VL

            # \\ Zustrommodell durchführen
            Speicherzustand = __Modell_Zustrom("oben", aktuellSekunden, ausgabezeit,
                                               Begleitdaten, dt_sub, Ausgabewerte,
                                               Speicherzustand)                         # neuen Speicherzustand nach Hinzufügen von Zellen oben berechnen

            # \\ Inversionen auflösen
            counterInv = 0
            Vp_zu = m_VL / __Modell_Stoffwerte("rho", T_Zustrom)
            start_inversion_status = __Modell_Inversionspruefung(Speicherzustand)       # Prüfen, ob Inversionen vorhanden sind
            inversion_status = start_inversion_status
            Speicherzustand = __Modell_Aufraumen(Speicherzustand)
            while inversion_status != "keine":                                          # Schleife läuft so lange, bis keine Inversionen mehr vorhanden sind
                counterInv += 1


                Speicherzustand = __Modell_Inversion(start_inversion_status,
                                                     "oben", Vp_zu, dt_sub,
                                                     Speicherzustand)                   # Inversionen auflösen
                Speicherzustand = __Modell_Aufraumen(Speicherzustand)                   # Speicher aufräumen
                if counterInv == 1:
                    Speicherzustand = __Modell_Impuls(start_inversion_status,
                                                       dt_sub, Speicherzustand)    # bei erstem Durchlauf der Schleife wird Impulsmodell durcgheführt
                    Speicherzustand = __Modell_Aufraumen(Speicherzustand)               # Speicher aufräumen
                inversion_status = __Modell_Inversionspruefung(Speicherzustand)         # Inversionsprüfung

            # \\ Horizontalmischung durchführen
            Speicherzustand = __Modell_Horizontalmischung(Speicherzustand)
            for hPos in list(Speicherzustand):                                          # Impuls auf Null setzen
                Speicherzustand[hPos][2] = 0
                Speicherzustand[hPos][3] = 0

            # \\ Masse und Energie im Speicher neu berechnen
            masse_zu = m_VL * dt_sub
            Masse_Global_Modellgrenzen += masse_zu
            Energie_Global_Modellgrenzen += masse_zu\
                                            * __Modell_Stoffwerte("h", T_Zustrom)

        # // Fall: Entladung durchführen - Abstrom

        elif m_VL < 0: # Massestrom im Vorlauf ist negativ
            Begleitdaten = {}
            Begleitdaten[aktuellSekunden] = {}
            Begleitdaten[aktuellSekunden]["m_Punkt"] = -m_VL
            Begleitdaten[aktuellSekunden]["T_Abstrom"] = T_Abstrom

            # \\ Abstrommodell durchführen
            Speicherzustand = __Modell_Abstrom("oben", aktuellSekunden,
                                               ausgabezeit, Ausgabewerte,
                                               Begleitdaten,dt_sub,
                                               Speicherzustand)                         # neuen Speicherzustand nach Verkleinern von Zellen oben berechnen

            # \\ Masse und Energie im Speicher neu berechnen
            masse_ab = m_VL * dt_sub
            theta_ab = Ausgabewerte["T_VL"][ausgabezeit]
            #theta_ab = T_Abstrom
            h_ab = __Modell_Stoffwerte("h", theta_ab)
            Masse_Global_Modellgrenzen += masse_ab
            Energie_Global_Modellgrenzen += masse_ab * h_ab
            H_Abstrom += -masse_ab * h_ab

            #abgefuerhte energie
            micro_energie = __Modell_Stoffwerte("cp", (T_Zustrom+theta_ab)/2 ) * m_VL * (T_Zustrom - theta_ab)
            macro_energie += micro_energie
        # ende if m_VL < 0
        # TODO counterInv ist eigenlich nicht noetig

        # // Zellgrößen anpassen (durch Teilen oder Zusammenlegen von Zellen)
        Speicherzustand = __Modell_Zellgroesse(j, "beides", Speicherzustand)

        # // Inversionen auflösen
        counterInv = 1
        start_inversion_status = __Modell_Inversionspruefung(Speicherzustand)           # auf Inversionen prüfen
        inversion_status = start_inversion_status
        while inversion_status != "keine":                                              # Schleife solange, bis alle Inversionen aufgelöst sind
            counterInv += 1
            Speicherzustand = __Modell_Inversion(start_inversion_status,
                                                 "nicht_definiert", 0, dt_sub,
                                                 Speicherzustand)                       # Inversionen auflösen
            Speicherzustand = __Modell_Aufraumen(Speicherzustand)                       # Speicher aufräumen
            inversion_status = __Modell_Inversionspruefung(Speicherzustand)             # auf Inversionen prüfen
        # ende while

        # // Waermeleitung berechnen
        # ACHTUNG Verschiedene Groessen vielleicht als input ermoeglichen
        # input koennte sein eine Config datei die immer da ist, und mit 
        # Standartwerte ausgefuellt ist
        # -------------
        theta_DR = speicher_param["T_DR"]                                               # Dampfraumtemperatur
        theta_o_Rand = Speicherzustand[max(list(Speicherzustand))][0]                   # Temperatur der obersten Schicht des Wassers im Speicher vor Ausführung der Wärmeleitung

        #Speicherzustand[max(list(Speicherzustand))][0] = theta_DR

        Fundamentzustand, Speicherzustand = __Modell_Waermeleitung(dt_sub,
                                                        theta_DR,
                                                        speicher_param["q_Punkt_U"],
                                                        Fundamentzustand,
                                                        Speicherzustand)                    # neuen Speicher- und Fundamentzustand nach Wärmeleitung berechnen
        all_h_pos = sorted(list(Speicherzustand))
        q_punkt_oben = (theta_DR - (theta_o_Rand+Speicherzustand[all_h_pos[-1]][0]) / 2)\
            / (Speicherzustand[all_h_pos[-1]][1]/2)\
            * __Modell_Stoffwerte("lambda",Speicherzustand[all_h_pos[-1]][0])               # Mittelwert Temp. vor und nach der Ausführung der Wärmeleitung, Länge ist halbe Zellhöhe der obersten Zelle, ergibt insgesamt den Wärmestrom bezogen auf die Querschnittsfläche des Speichers vom DR in die oberste Zelle
        # TODO -q_punkt_u + q_punkt_oben ?
        Energie_Global_Modellgrenzen += (-speicher_param["q_Punkt_U"] + q_punkt_oben)\
                                         * speicher_param["A_Quer"] * dt_sub                # neue Gesamtenergie des Speichers innerhalb der Modellgrenzen berechnen, indem Wärmezu/abfuhr im DR und Fundament? addiert wird
        BilanzEnergie -= (-speicher_param["q_Punkt_U"] + q_punkt_oben)\
                          * speicher_param["A_Quer"] * dt_sub                               # Frage: Was wird hier berechnet? Ist das eine Größe für die Änderung der Energie?
        Q_oben += q_punkt_oben * speicher_param["A_Quer"] * dt_sub                          # Wärmemenge, die vom DR an das Wasser abgegeben wird

        E_Verlust_Mantel, Kapazitaeten, Speicherzustand = __Modell_Kapazitaeten(
                                                                    dt_sub,
                                                                    T_amb,
                                                                    Kapazitaeten,
                                                                    Speicherzustand)        # Verluste über den Mantel berechnen (Wärmeübertragung zwischen Wasser und Mantel sowie Mantel und Umgebung)
        
        Speicherzustand = __Modell_Aufraumen(Speicherzustand)                               # kleine Zellen löschen, neue hPos berechnen
        E_Verlust_Mantel_alle_dt_sub += E_Verlust_Mantel                                    # Verlust über Mantel über alle Subzeitschritte aufaddieren
        Energie_Global_Modellgrenzen -= E_Verlust_Mantel                                    # Energieverlust innerhalb der Modellgrenzen berechnen
        # start_inversion_status = __Modell_Inversionspruefung(Speicherzustand)
        # inversion_status = start_inversion_status
        # while inversion_status != "keine":
        #     Speicherzustand = __Modell_Inversion(start_inversion_status, "nicht definiert",
        #                                          0, dt_sub, Speicherzustand)
        #     Speicherzustand = __Modell_Aufraumen(Speicherzustand)
        #     inversion_status = __Modell_Inversionspruefung(Speicherzustand)

        # for v in Speicherzustand.values():
        #     v[2] = 0
        #     v[3] = 0

        last = 1 if j==n_sub else 0                                                         # last wird 1, wenn letzter Zeitschritt erreicht ist, sonst 0

        # // Fuellstand und Bodendruck bestimmen
        # ------------------------------------
        all_h_pos = sorted(list(Speicherzustand))

        Ausgabewerte["FuellStand"][ausgabezeit] = all_h_pos[-1]\
                                            + Speicherzustand[all_h_pos[-1]][1] / 2         # Berechnung des Füllstands mithilfe der obersten Zelle des Speichers
        Ausgabewerte["p_unten"][ausgabezeit] = 0
        for hPos in all_h_pos:
            Ausgabewerte["p_unten"][ausgabezeit] += (g
                                                    * Speicherzustand[hPos][1] 
                            * __Modell_Stoffwerte("rho", Speicherzustand[hPos][0]))         # berechnet den Bodendruck in Pa, indem der Druck aller Zellen aufaddiert wird
        Ausgabewerte["p_unten"][ausgabezeit] /= 100000                                      # berechnet Bodendruck in bar
        
        masse_Speicher = __masse_berechnen(Speicherzustand)                                 # berechnet die Speichermasse in kg?
        Masse_Bilanz_Speicher = masse_Speicher - Masse_Global_Speicherzustand               # berechnet die Änderung der Speichermasse im Vergleich zur Masse zu Beginn des Zeitschritts
        #print(Masse_Bilanz_Speicher)
        #print(Masse_Global_Modellgrenzen)
        Masse_Bilanz_Korrektur = -Masse_Bilanz_Speicher + Masse_Global_Modellgrenzen        # Frage: unnötig, weil nicht verwendet?
        #print(Masse_Bilanz_Korrektur)

        # Das hier im Vollstaendigen Modell nicht vorhanden.
        # aber notwendig um die gesamtmasse konstant zu halten
        Speicherzustand[all_h_pos[0]][1] += ( Masse_Bilanz_Korrektur
                / (2 * speicher_param["A_Quer"]
                    * __Modell_Stoffwerte("rho", Speicherzustand[all_h_pos[0]][0])) )
        Speicherzustand[all_h_pos[-1]][1] += ( Masse_Bilanz_Korrektur
            / (2 * speicher_param["A_Quer"]
                * __Modell_Stoffwerte("rho", Speicherzustand[all_h_pos[-1]][0])) )


        # masse erneute berechnen da 2 zellen geaendert wurden
        masse_Speicher = __masse_berechnen(Speicherzustand)                                 # Frage: eigentlich nicht notwendig?
        Masse_Bilanz_Speicher = masse_Speicher - Masse_Global_Speicherzustand               # Frage: auch nicht notwendig?
        
        alle_H = [v[1] * __Modell_Stoffwerte("rho", v[0])
                       * __Modell_Stoffwerte("h",v[0]) 
                for k,v in sorted(Speicherzustand.items())]                                 # Enthalpie für jede Zelle berechnen
        alle_C_fundament = [v[0] * v[1]
                        * __Modell_Stoffwerte("rho_Fundament") 
                        * __Modell_Stoffwerte("cp_Fundament")  
                        for k,v in sorted(Fundamentzustand.items())]                        # Kapazitäten für jede Fundamentzelle berechnen
        alle_C_speicher = [v[0] * v[2] for k,v in sorted(Kapazitaeten.items())]             # Kapazitäten im Speichermantel berechnen

        # // Berechnung der Ausgabewerte
        Energie_Speicher = speicher_param["A_Quer"] * sum(alle_H)                           # Energie im Speicherwasser aus Enthalpien berechnen
        Energie_Speicher += speicher_param["A_Quer"] * sum(alle_C_fundament)                # Energie im Fundament addieren
        Energie_Speicher += sum(alle_C_speicher)
        e_ges = sum(alle_H) * speicher_param["A_Quer"]/1E09                                      # Energie im Speichermantel addieren
        Energie_Bilanz_Speicher = Energie_Speicher - Energie_Global_Speicherzustand         # Differenz zu Energie vor dem Zeitschritt berechnen
        Ausgabewerte["m_Abw_global_in_kg"][ausgabezeit] = Masse_Bilanz_Speicher\
                                                        - Masse_Global_Modellgrenzen        # Prüfen, ob im Speicher so viel Masse hinzugefügt wurde, wie durch Diff. m_zu und m_ab hätte hinzugefügt werden müssen
        Ausgabewerte["H_Abw_global_in_kJ"][ausgabezeit] = (Energie_Bilanz_Speicher
                                            - Energie_Global_Modellgrenzen) / 1000          # Prüfen, ob im Speicher so viel Energie hinzugefügt wurde, wie
        Ausgabewerte["check"][ausgabezeit] = speicher_param["A_Quer"] * sum(alle_H)         # Frage: Was wird hier gecheckt?
        #Speicherzustand = zone_finden(Speicherzustand)
        alle_Temperaturprofile = __Modell_Ausgabe_Zeitschritt(last,
                                                        ausgabezeit,
                                                        alle_Temperaturprofile,
                                                        Fundamentzustand,
                                                        Speicherzustand,
                                                        Kapazitaeten)                       # Temperaturprofile


    # // Rueckgabewerte bestimmen
    #----------------------------
    m_Abstrom = 0
    if m_VL < 0:
        m_Abstrom = -m_VL * dt
    if m_RL < 0:
        m_Abstrom = -m_RL * dt
    T_Abstrom = -1
    if m_Abstrom > 0:
        T_Abstrom = __Modell_Stoffwerte("h_rev", H_Abstrom / m_Abstrom)

    
    m_ges = masse_Speicher / 1000 # ausgabe in t
    E_ges  = Energie_Speicher / 1E09 # ausgabe in GJ
    # nutzbare masse
    # definiert als masse des Mediums die waermer als die definierte Grentztempe-
    # ratur ist und vollstaendig unterhalb der Oberkante des oberen Diffusors
    all_h_pos = sorted(list(Speicherzustand))
    lastKey = all_h_pos[-1]
    h_WS = lastKey + Speicherzustand[lastKey][1] / 2
    logger.debug("h_WS=%s", h_WS)
    h_pos_OK_oberer_Diff = h_WS - speicher_param["H_WS_OK_Dif"]
    logger.debug("h_pos_OK_oberer_Diff=%s", h_pos_OK_oberer_Diff)
    T_grenz = speicher_param["T_grenz"]
    m_nutz = __masse_nutz(Speicherzustand, h_WS)
    m_nutz_momentan = __masse_nutz(Speicherzustand, h_WS, T_RL)
    E_nutz = __energie_nutz(Speicherzustand, h_WS)
    E_nutz_momentan = __energie_nutz(Speicherzustand, h_WS, T_RL)
    if abs(m_RL) != 0:
        t_bis_leer = m_nutz*1000 / abs(m_RL)
    else:
        t_bis_leer = "inf"
    T_Diff_U = __Modell_Temperatur_Diffusorhoehe("unten", h_WS, Speicherzustand)
    T_Diff_O = __Modell_Temperatur_Diffusorhoehe("oben", h_WS, Speicherzustand)

    if T_Diff_O < T_grenz:
        warnings.warn("Die mittlere Temperatur am oberen Diffusor T_Diff_O ist kleiner als die Grenztemperatur")

    m_nutz_max = __masse_nutz_max(Speicherzustand, h_WS)

    mp_max_P_BL = speicher_param["Vp_max"] / 3600 * __Modell_Stoffwerte("rho", T_Diff_O)
    mp_max_P_EL = speicher_param["Vp_max"] / 3600 * __Modell_Stoffwerte("rho", T_Diff_U)

    mp_max_BL = min((m_nutz_max - m_nutz) * 1000 / dt , mp_max_P_BL)
    mp_max_BL = max(mp_max_BL, 0)
    mp_max_EL = min(m_nutz * 1000 / dt, mp_max_P_EL)

    mp_min_BL = speicher_param["Vp_min_rel"] * mp_max_P_BL
    mp_min_EL = speicher_param["Vp_min_rel"] * mp_max_P_EL

    mp_min = max(mp_min_EL, mp_min_BL)

    # // Verluste berechnen
    Q_V_DR = -Q_oben / dt # Q_oben wird berechnet mit Bezug auf Dampftemperatur (T_DR - T_medium) also DR erwaermt das Medium
    Q_V_Zyl = E_Verlust_Mantel_alle_dt_sub / dt
    Q_V_Erd = speicher_param["q_Punkt_U"] * speicher_param["A_Quer"]
    Q_V_ges = Q_V_Zyl + Q_V_Erd + Q_V_DR # Positive VERLUSTE

    # locating the thermocline 
    low_t = high_t = partial_volume_thermocline = charge_factor = thermocline_heights = None
    thermocline_min_h = None
    thermocline_max_h = None
    # // Rückgabewerte in eine Datei schreiben
    outputs = {
    "t" : t,
    "T_Austritt" : T_Abstrom,
    "m_nutz" : m_nutz,
    "m_nutz_momentan" : m_nutz_momentan,
    "m_nutz_max" : m_nutz_max,
    "E_nutz" : E_nutz,
    "E_nutz_momentan" : E_nutz_momentan,
    "m_ges" : m_ges,
    "enthalpie_alle" : e_ges,
    "E_ges" : E_ges,
    "T_Diff_U" : T_Diff_U,
    "T_Diff_O" : T_Diff_O,
    "Q_V_ges" : Q_V_ges,
    "Q_V_DR" : Q_V_DR,
    "Q_V_Zyl" : Q_V_Zyl,
    "Q_V_Erd" : Q_V_Erd,
    "t_bis_leer" : t_bis_leer,
    "H_WS" : h_WS,
    "mp_max_BL" : mp_max_BL,
    "mp_max_EL" : mp_max_EL,
    "mp_min" :  mp_min,
    "untere_temperatur_mischzone" : low_t,
    "obere_temperatur_mischzone" : high_t,
    "untere_hoehe_mischzone" : thermocline_min_h,
    "obere_hoehe_mischzone" : thermocline_max_h,
    "mischzone_groesse_relativ" : partial_volume_thermocline,
    "beladefaktor_nach_mischzone" : charge_factor,
    "speicherzustand" : Speicherzustand
    }

    #outputs = Speicherzustand

    return outputs

# // Ab jetzt folgen die weiteren definierten Funktionen:
# // ____________________________________________________


def __energie_nutz(sz: dict, h_ws: float, T_bezug: float = None) -> float:
    if T_bezug is None:
        T_bezug = speicher_param["T_grenz"]
        
    h_pos_OK_oberer_Diff = h_ws - speicher_param["H_WS_OK_Dif"]
    energie = 0.0
    bezugsenthalpie = __Modell_Stoffwerte("h", T_bezug)

    for k, v in sz.items():
        rho = __Modell_Stoffwerte("rho", v[0])
        h_diff = __Modell_Stoffwerte("h", v[0]) - bezugsenthalpie

        if v[0] > T_bezug:
            if (k + v[1] / 2) <= h_pos_OK_oberer_Diff:
                term = v[1] * rho * h_diff
                energie += term
            elif (k - v[1] / 2) < h_pos_OK_oberer_Diff:
                term = (h_pos_OK_oberer_Diff - (k - v[1] / 2)) * rho * h_diff
                energie += term

    E_nutz = speicher_param["A_Quer"] * energie / 1E09  # GJ
    return E_nutz

def __masse_nutz(speicherzustand: dict, h_WS: float, T_bezug: float = None) -> float:
    """
    Berechnet die nutzbare Masse im Speicher in Tonnen (d. h. Zelltemperatur liegt über Grenztemperatur und Zelle liegt unter der Oberkante des oberen Diffusors)
    speicherzustand: Speicherzustand
    h_WS: Höhe des Wasserspiegels
     nutzbare Masse in Tonnen
    """
    if T_bezug is None:
        T_bezug = speicher_param["T_grenz"]
    h_pos_OK_oberer_Diff = h_WS - speicher_param["H_WS_OK_Dif"]         # Höhenposition der Oberkante des oberen Diffusors
    m_nutz = 0.0
    debug_list = []
    for k, v in speicherzustand.items():
        if (v[0] >= T_bezug and k + v[1]/2 <= h_pos_OK_oberer_Diff and k - v[1]/2 >= speicher_param["H_B_UK_Dif"] ):
            m_debug = v[1] * __Modell_Stoffwerte("rho", v[0])
            m_nutz += v[1] * __Modell_Stoffwerte("rho", v[0])
        elif v[0] >= T_bezug and k - v[1]/2 < h_pos_OK_oberer_Diff and k + v[1]/2 > h_pos_OK_oberer_Diff:
            m_debug = (h_pos_OK_oberer_Diff - (k - v[1]/2)) * __Modell_Stoffwerte("rho", v[0])
            m_nutz += (h_pos_OK_oberer_Diff - (k - v[1]/2)) * __Modell_Stoffwerte("rho", v[0])
        elif v[0] >= T_bezug and k + v[1]/2 > speicher_param["H_B_UK_Dif"] and k - v[1]/2 < speicher_param["H_B_UK_Dif"]:
            m_debug = ((k + v[1]/2) - speicher_param["H_B_UK_Dif"]) * __Modell_Stoffwerte("rho", v[0])
            m_nutz += ((k + v[1]/2) - speicher_param["H_B_UK_Dif"] ) * __Modell_Stoffwerte("rho", v[0])
        else:
            m_debug = 0

        debug_list.append(m_debug)

    m_nutz_gesamt = speicher_param["A_Quer"] * m_nutz / 1000         # nutzbare Masse in Tonnen
    logger.debug("m_nutz_gesamt=%s", m_nutz_gesamt)
    return m_nutz_gesamt
def __masse_nutz_max(speicherzustand, h_WS):
    """
    Berechnet die Masse im Speicher in Tonnen, die maximal beladen werden kann (d. h. die Masse zwischen den Diffusoren)
    speicherzustand: Speicherzustand
    h_WS: Höhe des Wasserspiegels
     maximal nutzbare Masse in Tonnen
    """
    h_pos_OK_oberer_Diff = h_WS - speicher_param["H_WS_OK_Dif"]         # Höhenposition der Oberkante des oberen Diffusors

    m_plug_nutz = 0

    for k,v in speicherzustand.items():
        if (k - v[1]/2 >= speicher_param["H_B_UK_Dif"] and k + v[1]/2 <= h_pos_OK_oberer_Diff):
            m_plug_nutz += v[1] * __Modell_Stoffwerte("rho", v[0])
        elif k - v[1]/2 < h_pos_OK_oberer_Diff and k + v[1]/2 > h_pos_OK_oberer_Diff:
            m_plug_nutz += (h_pos_OK_oberer_Diff - (k - v[1]/2)) * __Modell_Stoffwerte("rho", v[0])
        elif k + v[1]/2 > speicher_param["H_B_UK_Dif"] and k - v[1]/2 < speicher_param["H_B_UK_Dif"]:
            m_plug_nutz += ((k + v[1]/2) - speicher_param["H_B_UK_Dif"]) * __Modell_Stoffwerte("rho", v[0])

    m_nutz_max = speicher_param["A_Quer"] * m_plug_nutz / 1000         # nutzbare Masse in Tonnen
    return m_nutz_max

def __Modell_Temperaturabsenkung_rohr(T_in, mp, t_amb, method="VDI_2055"):
    lambda_iso = 0.0275 # W/m*K waermeleitkoeff
    lambda_rohr = 52.33 # W/m*K
    d_i = 0.5618 #m
    dicke_rohr = 0.56E-2 #m
    rho_rohr = 7870 # kg/m3
    cp_rohr = 460 # J/kgK
    dicke_iso = 0.5E-1 #m
    alpha_i = 20 #W/m2K
    alpha_a = 4 #W/m2k
    rho_medium = 1000
    cp_medium = 4186 #J/kgK
    L = 260 #m
    d_a_iso = d_i + 2*dicke_rohr + 2*dicke_iso
    d_i_iso = d_a_iso - 2*dicke_iso
    d_a_rohr = d_i + 2*dicke_rohr
    pi = 3.1457

    if method == "VDI_2055":
        #k_x = 1 / (1/pi/alpha_i/d_i + 1/2/pi/lambda_rohr*log(d_a_rohr/d_i) + 1/2/pi/lambda_iso*log(d_a_iso/d_i_iso)+1/pi/d_a_iso/alpha_a)

        lambda_gesamt = ( (1/(2*pi)) * ( (1/lambda_rohr) * log(d_a_rohr/d_i) + (1/lambda_iso) * log(d_a_iso/d_i_iso)))**-1
        k_r = ( 1 / (pi * d_i * alpha_i) + 1 / (lambda_gesamt) + 1 / (pi * d_a_iso * alpha_a))**-1
        epsilon = k_r / (mp* cp_medium)
        t_out =  -(-T_in + t_amb) * exp(-epsilon*L) + t_amb
    return t_out



# // Speichermasse berechnen
def __masse_berechnen(Speicher):
    """
    Berechnet die Speichermasse (gesamte Wassermasse im Speicher inklusive des Wassers über bzw. unter den Diffusoren)
    Speicher: Speicherzustand
     Speichermasse in kg
    """
    masse_pro_m2 = sum([v[1] * __Modell_Stoffwerte("rho", v[0])
                        for v in Speicher.values()])                # Masse pro m² in kg/m²
    return speicher_param["A_Quer"] * masse_pro_m2                  # Masse in kg

# // Funktion: Initialisierung des Speichers im ersten Zeitschritt
def __Modell_Initialisierung(H_UEB=None, BeladeFaktor=None, iniZeitstempel=0.00, H_UEB_pos=None,
                             alle_Temperaturprofile={}, spline_uebernehmen=False,zustand={}, debug=False):
    """
    Initialisierung des Speichers im allerersten Zeitschritt
    H_UEB: Höhe des Übergangsbereichs
    BeladeFaktor: Prozentangabe, wie viel des Speichers beladen ist
    iniZeitstempel: Zeitstempel der Initialisierung
    H_UEB_pos: Position des Beginns des Übergangsbereichs von unten gesehen
    alle_Temperaturprofile: ??? Frage: wird aktuell nicht verwendet
    spline_uebernehmen: Spline für Temperaturverteilung übernehmen
    zustand: vorgegebenes Temperaturprofil für Speicher
    debug: Frage: ??
     alle_Temperaturprofile, Speicherzustand, Fundamentzustand, Kapazitaeten
    """
    Speicherzustand = {}
    Fundamentzustand = {}
    Kapazitaeten = {}
    alle_Temperaturprofile = {}
    dH = 0.05                                                                           # Frage: wird nicht verwendet
    if H_UEB is None:
        H_UEB = speicher_param["H_UEB_start"]                                           # falls H_UEB in Funktion nicht vorgegeben wurde, wird es hier gemacht
    if BeladeFaktor is None:
        try:
            BeladeFaktor = speicher_param["Beladefaktor_start"]                         # falls BeladeFaktor in Funktion nicht vorgegeben wurde, wird es hier gemacht
        except:
            BeladeFaktor = None
    if BeladeFaktor is None:
        H_UEB_pos = speicher_param["H_UEB_pos"]                                         # falls im config kein BeladeFaktor angegeben wurde, wird H_UEB_pos verwendet
    # Startzustand erzeugen
    f_a = 2 * 0.9061938 / H_UEB                                                         # Frage: wo kommen diese Zahlen her?
    if spline_uebernehmen == False:                                                     # wenn kein Temperaturprofil vorgegeben ist, sondern eins erstellt werden muss
        if H_UEB_pos is not None:                                                       # wenn die Position des Übergangsbereichs vorgegeben ist
            T_diff = speicher_param["T_oben"] - speicher_param["T_unten"]               # Temperaturdifferenz zwischen heißer (oben) und kalter (unten) Schicht
            T_mittel = (speicher_param["T_oben"] + speicher_param["T_unten"]) / 2       # mittlere Temperatur zwischen heißer und kalter Schicht
            dh = 0.05                                                                   # Zellenabstand
            hpos = 0
            T_Grad = (speicher_param["q_Punkt_U"]
                    * speicher_param["H_B_UK_Dif"] 
                    / __Modell_Stoffwerte("lambda", speicher_param["T_unten"]))         # Temperaturgradient, q_Punkt_U ist Wärmestromdichte!
            counter = speicher_param["H_B_UK_Dif"] / dh                                 # Anzahl Zellen unter der Unterkante des unteren Diffusors
            # Das Temperaturprofil vom Boden bis UK Dif erzeugen
            while counter > 0:                                                          # zählt bis alle Zellen unterhalb der UK des Diffusors erzeugt sind
                hpos += dh * 0.5                                                        # Position der Zelle (Zellmitte)
                theta = speicher_param["T_unten"] - T_Grad * dh * (counter - 0.5)       # Temperatur der Zelle
                Speicherzustand[hpos] = [theta, dh, 0, 0]                               # Speicherzustand schreiben
                hpos += dh * 0.5                                                        # Oberkante der Zelle
                counter -= 1

            # Profil weiter bis zum Anfang der UEB erzeugen
            counter = 1
            hpos_fix = hpos
            while hpos_fix + dh * counter < H_UEB_pos:                                  # solange Position unter dem Beginn des Übergangsbereichs ist
                hpos += dh * 0.5                                                        # Position der Zelle
                Speicherzustand[hpos] = [speicher_param["T_unten"], dh, 0, 0]           # Speicherzustand schreiben, Zelle hat untere Temperatur
                hpos += dh * 0.5                                                        # Oberkante der Zelle
                counter += 1
            dh_rest = H_UEB_pos - hpos                                                  # Abstand der Oberkante der letzten erstellten Zelle zum Beginn des Übergangsbereichs
            hpos += dh_rest * 0.5                                                       # Mitte zwischen Oberkante letzte Zelle und Beginn Übergangsbereich
            Speicherzustand[hpos] = [speicher_param["T_unten"], dh_rest, 0, 0]          # für diese Zelle Speicherzustand schreiben
            hpos += dh_rest * 0.5                                                       # hpos entspricht jetzt Beginn des Übergangsbereichs
            
            # UEB Profil erzeugen im "kalten Bereich"
            # * 0.5 weil die Haelfte der Zellen im "warmen" Bereich
            # Wegen erf() wir muessen von den counter absteigend zaehlen  
            counter = H_UEB / dh * 0.5                                                      # Anzahl Zellen im kalten Bereich
            while counter > 0:
                hpos += dh * 0.5 
                theta = T_mittel - T_diff * 0.5 * special.erf(f_a * dh * (counter - 0.5))
                Speicherzustand [hpos] = [theta, dh, 0, 0]
                hpos += dh * 0.5
                counter -= 1
            
            # UEB Profil im "warmen" Bereich. Counter aufsteigend
            counter = 1
            while counter < H_UEB / dh * 0.5:
                theta = T_mittel + T_diff * 0.5 * special.erf(f_a * dh * (counter - 0.5))
                hpos += dh * 0.5
                Speicherzustand[hpos] = [theta, dh, 0, 0]
                hpos += dh * 0.5
                counter += 1

            # Profil weiter bis Oberkante des oberen Diff
            counter = 1
            hpos_fix = hpos
            while hpos_fix + dh * counter < speicher_param["H_WS_max"] - speicher_param["H_WS_OK_Dif"]:
                hpos += dh * 0.5
                Speicherzustand[hpos] = [speicher_param["T_oben"], dh, 0, 0]
                hpos += dh * 0.5
                counter += 1
            dh_rest = speicher_param["H_WS_max"] - speicher_param["H_WS_OK_Dif"] - hpos
            hpos += dh_rest * 0.5
            Speicherzustand[hpos] = [speicher_param["T_oben"], dh_rest, 0, 0]
            hpos += dh_rest * 0.5

            # Profil von OK ob. Dif. bis WS
            counter = 1
            T_Grad = (speicher_param["T_DR"] - speicher_param["T_oben"])\
                    / speicher_param["H_WS_OK_Dif"]
            hpos_fix = hpos
            while hpos_fix + dh * counter < speicher_param["H_WS_max"]:
                theta = speicher_param["T_oben"] + T_Grad * dh * (counter - 0.5)
                hpos += dh * 0.5
                Speicherzustand[hpos] = [theta, dh, 0, 0]
                hpos += dh * 0.5
                counter += 1
            dh_rest = speicher_param["H_WS_OK_Dif"] - hpos
            theta = speicher_param["T_oben"] + T_Grad * dh_rest * (counter - 0.5)
            hpos += dh_rest * 0.5
            Speicherzustand[hpos] = [theta, dh_rest, 0, 0]
            hpos += dh_rest * 0.5
            
            
        else:                                                                               # wenn der Beladefaktor vorgegeben ist
            rho_T_max = __Modell_Stoffwerte("rho", speicher_param["T_oben"])
            H_aktiv = (speicher_param["H_WS_max"] 
                    - speicher_param["H_B_UK_Dif"]
                    - speicher_param["H_WS_OK_Dif"])

            m_beladen = H_aktiv * speicher_param["A_Quer"] * rho_T_max
            m_heiss = m_beladen * BeladeFaktor
            m_kalt = m_beladen - m_heiss
            T_diff = speicher_param["T_oben"] - speicher_param["T_unten"]
            T_mittel = (speicher_param["T_oben"] + speicher_param["T_unten"]) / 2
            dH = 0.05

            # kalte Seite des aktiven Speicherbereichs
            kalteSeite = {}
            counter = 1
            while m_kalt > 0:
                
                #factor = special.erf(f_a * dH * (counter - 0.5))
                if debug==False:
                    theta = T_mittel - T_diff / 2 * special.erf(f_a * dH * (counter - 0.5))
                else:
                    theta = T_mittel - T_diff/2
                rho = __Modell_Stoffwerte("rho", theta)
                if rho * dH * speicher_param["A_Quer"] <= m_kalt:
                    m_kalt -= rho * dH * speicher_param["A_Quer"]
                    kalteSeite[counter] = [theta, dH, 0, 0]
                else:
                    dH_Rest = m_kalt / rho / speicher_param["A_Quer"]
                    kalteSeite[counter] = [theta, dH_Rest, 0, 0]
                    m_kalt = 0
                counter += 1
                # with open("theta.dat", "a") as f:
                #     f.write("%.8f;  %.8f\n"%(theta,factor))
            # heisse Seite des aktiven Speicherbereichs
            heisseSeite = {}
            counter = 1
            while m_heiss > 0:
                if debug==False:
                    theta = T_mittel + T_diff / 2 * special.erf(f_a * dH * (counter - 0.5))
                else:
                    theta = T_mittel + T_diff/2
                rho = __Modell_Stoffwerte("rho", theta)

                if rho * dH * speicher_param["A_Quer"] <= m_heiss:
                    m_heiss -= rho * dH * speicher_param["A_Quer"]
                    heisseSeite[counter] = [theta, dH, 0, 0]
                else:
                    dH_Rest = m_heiss / rho / speicher_param["A_Quer"]
                    heisseSeite[counter] = [theta, dH_Rest, 0, 0]
                    m_heiss = 0
                counter += 1
            
            # Speicherbereich unterhalb der Unterkante des unteren Diffusors
            unterUnterkante = {}
            # vertikaler Temperaturgradient unterhalb der Unterkante des unteren
            # Diffusors in K/m
            T_Grad = (speicher_param["q_Punkt_U"] 
                    * speicher_param["H_B_UK_Dif"] 
                    / __Modell_Stoffwerte("lambda", speicher_param["T_unten"]))

            counter = 1
            while True:
                if dH * counter < speicher_param["H_B_UK_Dif"]:
                    theta = speicher_param["T_unten"] - T_Grad * dH * (counter - 0.5)
                    theta = speicher_param["T_unten"]
                    unterUnterkante[counter] = [theta, dH, 0, 0]


                else:
                    dH_Rest = speicher_param["H_B_UK_Dif"] - dH * (counter - 1)
                    theta = speicher_param["T_unten"] - T_Grad * (dH * counter - dH_Rest / 2)
                    theta = speicher_param["T_unten"]
                    unterUnterkante[counter] = [theta, dH_Rest, 0 ,0]
                    break
                counter += 1
            
            # Speicherbereich oberhalb der Oberkante des oberen Diffusors
            ueberOberkante = {}
            # vertikaler Temperaturgradient oberhalb der Oberkante des oberen
            # Diffusors in K/m
            T_Grad = (speicher_param["T_DR"] - speicher_param["T_oben"])\
                    / speicher_param["H_WS_OK_Dif"]
            counter = 1
            while True:
                if dH * counter < speicher_param["H_WS_OK_Dif"]:
                    theta = speicher_param["T_oben"] + T_Grad * dH * (counter - 0.5)
                    ueberOberkante[counter] = [theta, dH, 0, 0]
                else:
                    dH_Rest = speicher_param["H_WS_OK_Dif"] - dH * (counter - 1)
                    theta = speicher_param["T_oben"] + T_Grad * (dH * counter - dH_Rest / 2)
                    ueberOberkante[counter] = [theta, dH_Rest, 0, 0]
                    break
                counter += 1
            # alles zusammensetzen
            h_Pos = 0
            for counter in sorted(list(unterUnterkante), reverse=True):
                h_Pos += unterUnterkante[counter][1] / 2
                Speicherzustand[h_Pos] = unterUnterkante[counter]
                h_Pos += unterUnterkante[counter][1] / 2

            for counter in sorted(list(kalteSeite), reverse=True):
                h_Pos += kalteSeite[counter][1] / 2
                Speicherzustand[h_Pos] = kalteSeite[counter]
                h_Pos += kalteSeite[counter][1] / 2

            for counter in sorted(list(heisseSeite)):
                h_Pos += heisseSeite[counter][1] / 2
                Speicherzustand[h_Pos] = heisseSeite[counter]
                h_Pos += heisseSeite[counter][1] / 2

            for counter in sorted(list(ueberOberkante)):
                h_Pos += ueberOberkante[counter][1] / 2
                Speicherzustand[h_Pos] = ueberOberkante[counter]
                h_Pos += ueberOberkante[counter][1] / 2

    else:
        # Fallunterscheidung wenn p nicht gegeben ist und einfach
        # interpoliert werden soll
        if speicher_param["p_unten_anfang"] != 0:                                                                                   # wird durchlaufen, wenn ein Startspeicherzustand vorgegeben wurde
            h_WS = speicher_param["p_unten_pos"]                                                # Höhe der Wassersäule wird auf Position des Fülldrucksensors gesetzt
            bgd_p = speicher_param["p_unten_anfang"] * 100000                                   # entspricht Anfangsfülldruck in Pa
            hPos_lst = sorted(list(zustand))                                                    # Liste, in der die Höhenpositionen aus "zustand" stehen, Frage: Wie sieht "zustand" aus?
            dh = hPos_lst[1] - hPos_lst[0]                                                      # Abstand der ersten beiden Höhenpositionen - MUSS Abstand zwischen allen Zellen sein!
            # Wasserspiegelhoehe bestimmen
            for i in range(0,len(hPos_lst)):
                
                if hPos_lst[i] < h_WS:                                                          # nächstes Schleifenelement, wenn aktuelle Höhenposition unter h_WS liegt, zu Anfang ist dies die Position des Fülldrucksensors
                    continue
                if bgd_p<0:                                                                     # nächstes Schleifenelement, wenn Fülldruck "aufgebraucht ist"
                    continue
                rho = __Modell_Stoffwerte("rho", zustand[hPos_lst[i]])                          # Dichte in Zelle berechnen
                bgd_p -= rho * g * (hPos_lst[i]-h_WS)                                           # Fülldruck um den Druck der "Zelle" reduzieren
                h_WS = hPos_lst[i]                                                              # h_WS wird auf aktuelle hPos gesetzt
                if bgd_p < 0:
                    h_WS += bgd_p / (g* rho)                                                    # falls der Fülldruck jetzt "aufgebraucht ist", wird die Wasserspiegelhöhe um das "zu viel" reduziert
            #h_WS = 38
            # Speicherzustand von unterste DTS bis Wasserspiegel befuehlen
            for i in range(0, len(hPos_lst)):
                if hPos_lst[i] + dh/2 <= h_WS:                                                  # für alle Zellen, die komplett unter dem Wasserspiegel liegen
                    Speicherzustand[hPos_lst[i]] = [zustand[hPos_lst[i]], dh, 0, 0]             # Speicherzustand wird für diesen Punkt erstellt
                    #print(hPos_lst[i])
                elif hPos_lst[i] - dh/2 <=h_WS:                                                 # für oberste Zelle, die teilweise unter dem Wasserspiegel liegt
                    dh_oberste = dh/2 - (hPos_lst[i] - h_WS)                                    # Höhe der Zelle, die noch unter dem Wasserspiegel liegt
                    hPos_oberste = hPos_lst[i] - dh/2 + dh_oberste/2                            # mittlere Position der neuen Zelle unter dem Wasserspiegel berechnen
                    Speicherzustand[hPos_oberste] = [zustand[hPos_lst[i]], dh_oberste, 0, 0]    # Speicherzustand wird für oberste Zelle geschrieben
                else:
                    break                                                                       # bei erster Zelle, die höher als Wasserspiegel liegt, wird abgebrochen
            #print(hPos_oberste)
            # Speicherzustand vom Boden bis unterste DTS befuehlen
            theta_grad = (Speicherzustand[min(hPos_lst)][0] - speicher_param["T_Boden"])\
                        / (min(hPos_lst) - speicher_param["H_Bodenstrecke"])                    # Temperaturgradient zwischen unterster DTS-Messstelle und Boden
            H_ueber_Boden = min(hPos_lst) - dh/2                                                # unterer Rand der untersten DTS-Messstellen-Zelle
            while H_ueber_Boden > 0:
                if H_ueber_Boden - dh >= 0:                                                     # für Zelle, die komplett über dem Boden liegt
                    hPos = H_ueber_Boden - dh/2                                                 # hPos der Zelle
                    theta = zustand[min(hPos_lst)] - theta_grad * (min(hPos_lst) - hPos)        # Temperatur der Zelle wird interpoliert
                    Speicherzustand[hPos] = [theta, dh, 0 , 0]                                  # Speicherzustand schreiben
                elif H_ueber_Boden > 0:                                                         # für letzte Zelle, die noch teilweise über dem Boden liegt
                    dh_unterste = H_ueber_Boden                                                 # Zellhöhe entspricht der Resthöhe
                    hPos = dh_unterste / 2                                                      # hPos ist Mitte der Zelle
                    theta = zustand[min(hPos_lst)] - theta_grad * (min(hPos_lst) - hPos)        # Temperatur wird interpoliert
                    Speicherzustand[hPos] = [theta, dh_unterste, 0, 0]                          # Speicherzustand schreiben
                else:
                    break                                                                       # Abbruch, sobald unter dem Boden
                H_ueber_Boden -= dh
                                                                     # unterer Rand der aktuellen Zelle für nächsten Schleifendurchlauf
        else:
            zustand[0.0] = speicher_param["T_Boden"]
            zustand[speicher_param["H_WS_max"]] = speicher_param["T_DR"]
            z = np.array(sorted(zustand.keys()))
            T = np.array([zustand[k] for k in z])
            z_new = np.arange(0.0, speicher_param["H_WS_max"], 0.1)
            interp = interpolate.interp1d(z, T, kind="linear")
            T_new = interp(z_new)
            dh = float(z_new[1] - z_new[0])
            Speicherzustand = {
                float(z) : [float(T), dh, 0, 0]
                for z, T in zip(z_new, T_new)
            }
    # for hPos in sorted(list(Speicherzustand)):
    #     with open("temp_init.dat", "a") as f:
    #         f.write("%.3f;  %.5f;  %.5f\n"%(hPos, Speicherzustand[hPos][1], Speicherzustand[hPos][0]))
    Speicherzustand = __Modell_Zellgroesse(0, "beides", Speicherzustand)                    # Zellen teilen bzw. zusammenlegen
    #// Fundamentzustand
    #-------------------
    all_h_pos  = sorted(list(Speicherzustand))
    theta_Grad_bodennah = (Speicherzustand[all_h_pos[1]][0] 
                           - Speicherzustand[all_h_pos[0]][0]) \
                           / (all_h_pos[1] - all_h_pos[0])                                  # Temperaturgradient zwischen den beiden untersten Speicherzellen
    theta_Kontakt = Speicherzustand[all_h_pos[0]][0] \
                    - Speicherzustand[all_h_pos[0]][1] / 2 * theta_Grad_bodennah            # extrapolierte Temperatur am Boden des Speichers
    theta_Erdreich = speicher_param["T_Erdreich"]                                                                     # Bodentemperatur
    q_punkt_Erdreich = speicher_param["q_Punkt_U"]                                                                   # Wärmestromdichte = Verlust an den Boden (Wärmestrom bezogen auf die Fläche) W/m²
    lambda_Fundament = speicher_param["lambda_fundament"]                                                                    # Wärmeleitfähigkeit Fundament W/(mK)
    h_Fundament = lambda_Fundament \
                  * (theta_Kontakt - theta_Erdreich) \
                  / q_punkt_Erdreich                                                        # Höhe des Fundment wird BERECHNET! Frage: ist das so gewollt? Kennen wir nicht die exakte Höhe?
    n_Fundament = 500                                                                       # Stützstellen im Fundament = Zellenanzahl
    dh_Fundament = h_Fundament / n_Fundament                                                # Zellenhöhe
    theta_Grad_Fundament = q_punkt_Erdreich / lambda_Fundament                              # Temperaturgradient im Fundament in K/m
    for i in range(0,n_Fundament):
        hPos = -(i + 0.5) * dh_Fundament                                                    # Positionen der Zellen
        theta = theta_Kontakt + hPos * theta_Grad_Fundament                                 # Temperatur der Zellen
        Fundamentzustand[hPos] = [theta, dh_Fundament]                                      # Fundamentzustand schreiben
    

    # // Zustand Kapazitaeten feste Einbauten und Wand
    # Kapazitaeten sind nicht spezifisch!
    # ---------------------------------------------
    n_steps = 1000                                                                          # Anzahl Zellen im Mantel
    dh_kapa = speicher_param["H_Mantel"] / n_steps                                          # Zellenhöhe im Mantel
    Massen_pro_m = pi *__Modell_Stoffwerte("rho_Mantel") \
                   * ((speicher_param["R_innen"] + 0.016)**2 
                      - speicher_param["R_innen"]**2)                                       # Masse pro Höhenmeter des Mantels
    kapazitaet = dh_kapa * Massen_pro_m * __Modell_Stoffwerte("cp_Mantel")                  # Wärmekapazität = Q / delta_T
    
    # Spline aus vertikalen Temperaturfeld bestimmen
    all_theta = [v[0] for k,v in sorted(Speicherzustand.items())]                           # Temperaturen aus Speicherzustand in eine Variable schreiben
    spline = interpolate.CubicSpline(all_h_pos,all_theta,bc_type="natural")                 # Spline für Temperatur im Speicher erstellen
    for i in range(0, n_steps):
        hPos = (i + 0.5) * dh_kapa                                                          # Position der Mantelzellen
        if hPos < all_h_pos[-1]:
            theta = float(spline(hPos))                                                            # wenn Mantelzelle unter Wasserspiegel liegt, wird Temperatur aus Spline verwendet
        else:
            theta = float(all_theta[-1])                                                           # über dem Wasserspiegel wird Temperatur auf Temperatur der obersten Speicherzelle gesetzt
        Kapazitaeten[hPos] = [theta, dh_kapa, kapazitaet]                                   # Kapazitäten schreiben

    alle_Temperaturprofile = __Modell_Ausgabe_Zeitschritt(1, iniZeitstempel,
                                                          alle_Temperaturprofile,
                                                          Fundamentzustand,
                                                          Speicherzustand,
                                                          Kapazitaeten)                     # Temperaturprofil für Initialzustand in alle_Temperaturprofile ablegen

    #return Speicherzustand
    return alle_Temperaturprofile, Speicherzustand, Fundamentzustand, Kapazitaeten

# // Funktion: Zellgrößen anpassen
def __Modell_Zellgroesse(step, teilenZusammen, Speicherzustand):
    """
    Große Zellen dritteln (Randzellen halbieren) und zu kleine Zellen mit der Zelle darüber (oberste Zelle darunter) zusammenlegen.
    step: Frage: was ist das?
    teilenZusammen: teilen: nur teilen, zusammen: nur zusammenlegen, beides: teilen und zusammenlegen
    Speicherzustand: Speicherzustand
     Speicherzustand
    """
    dh_max = speicher_param["max_cell_height"]                                                                    # maximal zulässige Zellenhöhe

    # maximal zullaessige zellenhoehe * Delta_theta

    if teilenZusammen == "teilen":
        dh_max_theta = 1.5E-01                                                      # Frage: was gibt das an?
    else:
        dh_max_theta = 1.5E-02                                                      # Frage: same here?
    if teilenZusammen in ["beides", "teilen"]:
        zuGrossezelle = True
        while zuGrossezelle:                                                        # solange durchführen, wie es zu große Zellen gibt -> wird erneut durchgeführt, wenn mindestens eine Zelle zu groß war
            Speicherzustand_neu = {}                                                # Variable für neuen Speicherzustand
            all_h_pos = sorted(list(Speicherzustand))
            zuGrossezelle = False
            for i in range(1, len(all_h_pos)-1):                                    # alle bis auf die oberste und unterste Zelle durchgehen
                hPos = all_h_pos[i]                                                 # Position der Zelle selbst
                hPos_minus = all_h_pos[i-1]                                         # Position der Zelle darunter
                hPos_plus = all_h_pos[i+1]                                          # Position der Zelle darüber
                theta_i = Speicherzustand[hPos][0]                                  # Temperatur der zelle selbst
                theta_minus = Speicherzustand[hPos_minus][0]                        # Temperatur der Zelle darunter
                theta_plus = Speicherzustand[hPos_plus][0]                          # Temperatur der Zelle darüber
                dh_i = Speicherzustand[hPos][1]                                     # Größe der Zelle selbst
                
                if ((dh_i > dh_max)                                                 # wenn Zelle größer als maximal erlaubte Höhe ist
                or ((abs(theta_i - theta_minus))**0.75 * dh_i > dh_max_theta)       # oder die Temperaturdifferenz zwischen dieser und der unteren Zelle zu groß ist
                or ((abs(theta_plus - theta_i))**0.75 * dh_i > dh_max_theta)):      # oder die Temperaturdifferenz zwischen dieser und der oberen Zelle zu groß ist
                    # with open("pyCheck.dat", "a") as f:
                    #     f.write("%.4f; %.8f; %d  |||  %.3f;  %.3f;  %.3f\n"%(len(Speicherzustand_neu), dh_i, i, theta_minus, theta_i, theta_plus))
                    dh_minus = (dh_i + Speicherzustand[hPos_minus][1]) / 2          # Abstand zwsch. Zellmittelpunkten mit Zelle darunter - entspricht der halben Zellhöhe beider Zellen zusammen
                    dh_plus = (dh_i + Speicherzustand[hPos_plus][1]) / 2            # Abstand zwsch. Zellmittelpunkten mit Zelle darüber - entspricht der halben Zellhöhe beider Zellen zusammen

                    grad_minus = (theta_i - theta_minus) / dh_minus                 # Temperaturgradient mit Zelle darunter
                    grad_plus = (theta_plus - theta_i) / dh_plus                    # Temperaturgradient mit Zelle darüber

                    theta_u = theta_i - grad_minus * dh_i / 3                       # Temperatur 1/3 der Zellhöhe unter der Zellmitte -> Zelle soll gedrittelt werden
                    theta_o = theta_i + grad_plus * dh_i / 3                        # Temperatur 1/3 der Zellhöhe über der Zellmitte

                    m_t = dh_i * __Modell_Stoffwerte("rho", theta_i)                # Masse der Zelle bezogen auf die Speicherfläche
                    h_t = __Modell_Stoffwerte("h", theta_i)                         # Enthalpie der Zelle
                    loop = 0
                    while loop < 2:                                                 # Loop zweimal durchführen
                        # kann arithmetisch gemittelt werden, da alle Massen gleich
                        d_h = h_t - (__Modell_Stoffwerte("h", theta_u) 
                                     + __Modell_Stoffwerte("h",theta_o) 
                                     + __Modell_Stoffwerte("h", theta_i)) / 3       # Abweichung der "neuen" Enthalpie von der aktuellen Enthalpie der Zelle
                        d_theta = d_h / __Modell_Stoffwerte("cp", theta_i)          # Temperaturdifferenz zum Ausgleich der Enthalpiedifferenz
                        theta_i += d_theta                                          # Korrektur der Temperatur
                        theta_u += d_theta                                          # Korrektur der Temperatur
                        theta_o += d_theta                                          # Korrektur der Temperatur
                        loop += 1

                    dh_u = m_t / 3 / __Modell_Stoffwerte("rho", theta_u)            # Berechnung der Zellhöhe der neuen unteren Zelle
                    dh_o = m_t / 3 / __Modell_Stoffwerte("rho", theta_o)            # Berechnung der Zellhöhe der neuen oberen Zelle
                    dh_i = m_t / 3 / __Modell_Stoffwerte("rho", theta_i)            # Berechnung der Zellhöhe der neuen mittleren Zelle
                    if teilenZusammen == "teilen":                                  # falls nur geteilt werden soll -> es sind Impulskomponenten vorhanden
                        Speicherzustand_neu[hPos - dh_i / 3] = [
                                                    theta_u,
                                                    dh_u,
                                                    Speicherzustand[hPos][2],
                                                    Speicherzustand[hPos][3]
                                                    ]                               # neuen Speicherzustand für untere Zelle schreiben

                        Speicherzustand_neu[hPos + dh_i / 3] = [
                                                    theta_o,
                                                    dh_o,
                                                    Speicherzustand[hPos][2],
                                                    Speicherzustand[hPos][3]
                                                    ]                               # neuen Speicherzustand für obere Zelle schreiben

                        Speicherzustand_neu[hPos] = [
                                                    theta_i,
                                                    dh_i,
                                                    Speicherzustand[hPos][2],
                                                    Speicherzustand[hPos][3]
                                                    ]                               # neuen Speicherzustand der mittleren Zelle schreiben
                    else:                                                           # falls geteilt und zusammengelegt werden soll -> Impulskomponenten sind Null
                        Speicherzustand_neu[hPos - dh_i / 3] = [theta_u, dh_u, 0, 0] # neuen Speicherzustand für untere Zelle schreiben
                        Speicherzustand_neu[hPos + dh_i / 3] = [theta_o, dh_o, 0, 0] # neuen Speicherzustand für obere Zelle schreiben
                        Speicherzustand_neu[hPos] = [theta_i, dh_i, 0, 0]            # neuen Speicherzustand für mittlere Zelle schreiben
                    zuGrossezelle = True

                    # Ende if
                else:                                                               # Zelle ist nicht zu groß
                    Speicherzustand_neu[hPos] = Speicherzustand[hPos]               # neuer Speicherzustand der Zelle entspricht altem Zustand
                    # Ende else
                # Ende for
            # Randzellen
            for i in [0, len(all_h_pos)-1]:                                         # oberste und unterste Speicherzelle
                hPos = all_h_pos[i]
                dh_i = Speicherzustand[hPos][1]
                if dh_i > dh_max:
                    Speicherzustand_neu[hPos - dh_i / 4] = [
                                                    Speicherzustand[hPos][0],
                                                    dh_i / 2,
                                                    0,
                                                    0]                              # Zelle wird halbiert und beide Zellen erhalten die Temperatur der alten Zelle

                    Speicherzustand_neu[hPos + dh_i / 4] = [
                                                    Speicherzustand[hPos][0],
                                                    dh_i / 2,
                                                    0,
                                                    0]
                    zuGrossezelle = True

                else:
                    Speicherzustand_neu[hPos] = Speicherzustand[hPos]
                # Ende for
            Speicherzustand = __Modell_Aufraumen(Speicherzustand_neu)               # kleine Zellen löschen und Positionen neu bestimmen
        # ende while

    #ende If

    # zu kleine Zellen zusammenlegen
    if teilenZusammen in ["beides", "zusammen"]:
        f_hyst = 4                                                                  # dient dazu, dass Zellen, die zuvor geteilt wurden, nicht wieder vereint werden
        if step < 1:
            f_hyst = 2
        all_h_pos = sorted(list(Speicherzustand))
        count = 0

        for i in range(0, len(Speicherzustand)):                                    # alle Zellen durchgehen
            # an dieser Stelle hPos_i == hPos also hPos lassen (vgl. Perl)
            hPos = all_h_pos[i]
            theta_i = Speicherzustand[hPos][0]
            dh_i = Speicherzustand[hPos][1]
            mass_i = dh_i * __Modell_Stoffwerte("rho", theta_i)                     # Masse bezogen auf die Speicherfläche

            if i < len(all_h_pos)-1:                                                # für alle Zellen außer die oberste
                hPos_dazu = all_h_pos[i+1]                                          # Zelle darüber wird "addiert"
            if i == len(Speicherzustand)-1:                                         # für die oberste Zelle
                count_mix = 0
                while Speicherzustand[all_h_pos[i-1-count_mix]][1] == 0:            # Finden der ersten Zelle unter der obersten Zelle, die eine Höhe größer als Null hat
                    count_mix += 1
                    if i-1-count < 0:                                               # das passiert, falls jede vorherige Zelle mit der Zelle darüber zusammengelegt wurde, dann ist i = count
                        logger.warning("Keine Zelle zum Mischen gefunden. Abbruch!")
                        break
                hPos_dazu = all_h_pos[i - 1 - count_mix]                            # Position der Zelle, die "addiert" wird
            theta_dazu = Speicherzustand[hPos_dazu][0]                              # Temperatur der Zelle, die "addiert" werden soll
            # die beiden Zahlenfaktoren muessen hinreichen klein sein, damit zuvor
            # geteilte Zellen (s. weiter unter) nicht wieder zusammen gelegt werden
            if ((abs(theta_dazu - theta_i)**0.75 * dh_i < dh_max_theta / f_hyst)    # Temperaturgradient muss hinreichend klein sein
            and (dh_i < dh_max / f_hyst)):                                          # und Zelle kleiner als eine bestimmte Höhe sein
                mass_dazu = (Speicherzustand[hPos_dazu][1]
                             * __Modell_Stoffwerte("rho", theta_dazu))              # Masse der Zelle, die addiert wird

                H_neu = (mass_i * __Modell_Stoffwerte("h", theta_i)
                         + mass_dazu * __Modell_Stoffwerte("h", theta_dazu))        # neue Enthalpie der Zelle

                theta_neu = __Modell_Stoffwerte("h_rev", H_neu / (mass_i+mass_dazu))    # neue Temperatur der Zelle
                dh_neu = (mass_i+mass_dazu) /  __Modell_Stoffwerte("rho", theta_neu)    # neue Höhe der Zelle
                Speicherzustand[hPos] = [0, 0, 0, 0]                                    # Speicherzustand der betrachteten Zelle auf Null setzen, sodass diese beim Aufräumen gelöscht wird
                Speicherzustand[hPos_dazu] = [theta_neu, dh_neu, 0, 0]                  # Speicherzustand der addierten Zelle neu setzen
                count += 1

    Speicherzustand = __Modell_Aufraumen(Speicherzustand)                               # kleine Zellen löschen und Positionen der Zellen neu bestimmen

    return Speicherzustand

# // Funktion: Aufräumen
def __Modell_Aufraumen(Speicherzustand):
    """
    Entfernen von sehr kleinen Zellen sowie Neubestimmung von hPos
    Speicherzustand: Speicherzustand
     Speicherzustand
    """

    hPosNeu = 0
    SpeicherzustandAlt = Speicherzustand                                                # aktueller Speicherzustand wird in diese Variable geschrieben
    Speicherzustand = {}                                                                # Speicherzustand wird gelöscht
    for hPosAlt in sorted(list(SpeicherzustandAlt)):
        if SpeicherzustandAlt[hPosAlt][1] < 1.0E-09:
            continue                                                                    # sehr kleine Zellen werden gelöscht
        else:

            hPosNeu += SpeicherzustandAlt[hPosAlt][1] / 2                               # zu hPosNeu wird von Null beginnend die halbe Zellhöhe addiert = entspricht der Zellmitte
            Speicherzustand[hPosNeu] = SpeicherzustandAlt[hPosAlt]                      # alter Speicherzustand wird neuer hPos zugewiesen
            hPosNeu += SpeicherzustandAlt[hPosAlt][1] / 2                               # erneut halbe Zellhöhe addiert = entspricht der oberen Zellgrenze

    return Speicherzustand

# // Modell: Horizontalmischung
def __Modell_Horizontalmischung(Speicherzustand):
    all_h_pos = sorted(list(Speicherzustand))

    for i in range(0,len(all_h_pos)):
        hPos = all_h_pos[i]
        if Speicherzustand[hPos][3] <=0 :
            continue
        counter_plus = 1
        counter_minus = 1
        mix_plus = 1
        mix_minus = 1
        v_wirk = Speicherzustand[hPos][3]
        dh_drho_plus = 0
        dh_drho_minus = 0
        dh_plus_sum = 0
        dh_minus_sum = 0
        dh_max = tan(6.5/180*pi)\
                 * (speicher_param["R_innen"] - speicher_param["R_Dif"]) 
        
        oberste_start = 0
        unterste_start = 0

        if i >= len(all_h_pos) - 2:
            oberste_start = 1
            counter_plus = 0
        if i <= 1:
            unterste_start = 1
            counter_minus = 0

        while mix_plus==1 or mix_minus==1:
            hPos_plus = all_h_pos[i+counter_plus]
            hPos_minus = all_h_pos[i-counter_minus]
            dh_plus = Speicherzustand[hPos_plus][1]
            dh_plus_sum += dh_plus
            dh_minus = Speicherzustand[hPos_minus][1]
            dh_minus_sum += dh_minus
            theta_plus = Speicherzustand[hPos_plus][0]
            theta_minus = Speicherzustand[hPos_minus][0]

            rho = __Modell_Stoffwerte("rho", Speicherzustand[hPos][0])
            rho_plus = __Modell_Stoffwerte("rho", theta_plus)
            rho_minus = __Modell_Stoffwerte("rho", theta_minus)

            dh_pot_plus = 0
            # B: ACHTUNG: bitte checken ob das eigentlich so ist (unless)
            if rho == rho_plus:
                if not(oberste_start==1 or mix_plus==0):
                    dh_pot_plus = dh_plus
            else:
                if not (oberste_start==1 or mix_plus==0):
                    dh_pot_plus = (rho / (2 * g) *v_wirk**2 - dh_drho_plus)\
                                /abs(rho - rho_plus)
            dh_pot_minus = 0
            if rho == rho_minus:
                if not (unterste_start==1 or mix_minus==0):
                    dh_pot_minus = dh_minus
            else:
                if not (unterste_start==1 or mix_minus==0):
                    dh_pot_minus = (rho / (2 * g) *v_wirk**2 - dh_drho_minus)\
                                    /abs(rho - rho_minus)
            dh_pot_minus = min(dh_pot_minus, dh_minus)
            dh_pot_plus = min(dh_pot_plus, dh_plus)

            mix_plus = 0
            mix_minus = 0
            if dh_pot_plus > 1.0E-09:
                mix_plus = 1
            if dh_pot_minus > 1.0E-09:
                mix_minus= 1
            
            # B: ACHTUNG: abweichung von perl check if richtig
            #FIXED
            if dh_plus_sum > dh_max:
                mix_plus = 0
            if dh_minus_sum > dh_max:
                mix_minus = 0
            if oberste_start == 1:
                mix_plus = 0
            if unterste_start == 1:
                mix_minus = 0
            if (i + counter_plus == len(all_h_pos) - 2) and (dh_plus == 0):
                mix_plus = 0
            if (i - counter_minus == 1) and (dh_minus == 0):
                mix_minus = 0

            dh_drho_plus += 1 * dh_pot_plus * abs(rho - rho_plus)
            # B: ACHTUNG: Abweichung vom PERL
            f_minus = (hPos_minus / (1.5 * (speicher_param["H_B_UK_Dif"]+speicher_param["H_RS_Dif"]))) ** 2
            f_minus = min(f_minus, 1)
            f_minus = max(f_minus, 0.1)
            dh_drho_minus += f_minus * dh_pot_minus * abs(rho - rho_minus)

            if mix_plus == 1:
                dh = Speicherzustand[hPos][1]
                theta = Speicherzustand[hPos][0]
                rho = __Modell_Stoffwerte("rho", theta)
                m = dh * rho
                m_plus = rho_plus * dh_pot_plus
                H = m * __Modell_Stoffwerte("h", theta)\
                    + m_plus * __Modell_Stoffwerte("h", theta_plus)
                Speicherzustand[hPos][0] = __Modell_Stoffwerte("h_rev", (H/(m+m_plus)) )
                rho = __Modell_Stoffwerte("rho",Speicherzustand[hPos][0])
                Speicherzustand[hPos][1] = (m + m_plus) / rho
                Speicherzustand[hPos][3] = Speicherzustand[hPos][3]*m / (m+m_plus)
                Speicherzustand[hPos_plus][1] -= dh_pot_plus
                if i + counter_plus < len(all_h_pos) - 2:
                    counter_plus += 1
            
            if mix_minus == 1:
                dh = Speicherzustand[hPos][1]
                theta = Speicherzustand[hPos][0]
                rho = __Modell_Stoffwerte("rho", theta)
                m = dh * rho
                m_minus = rho_minus * dh_pot_minus
                H = m * __Modell_Stoffwerte("h", theta)\
                    + m_minus * __Modell_Stoffwerte("h", theta_minus)
                Speicherzustand[hPos][0] = __Modell_Stoffwerte("h_rev",
                                                               H / (m + m_minus))
                rho = __Modell_Stoffwerte("rho",Speicherzustand[hPos][0])
                Speicherzustand[hPos][1] = (m + m_minus) / rho
                Speicherzustand[hPos][3] = Speicherzustand[hPos][3]\
                                           * m / (m + m_minus)
                Speicherzustand[hPos_minus][1] -= dh_pot_minus
                if i - counter_minus > 1:
                    counter_minus += 1     
       # ende While
        break
    # end for i in range()
    Speicherzustand = __Modell_Aufraumen(Speicherzustand)
    return Speicherzustand

# // Funktion: Impuls
def __Modell_Impuls(inversion_status, Zeitabstand, Speicherzustand):

    all_h_pos = sorted(list(Speicherzustand))
    # empirisch
    f_imp_an_b =  Zeitabstand / 3

    for k,v in Speicherzustand.items():
        if v[2] < 4E-03:
            v[2] = 0
    indices = [i for i in range(0, (len(Speicherzustand) - 1) )]
    if inversion_status == "fallend":
        Speicherzustand[all_h_pos[0]][2] = 0
    elif inversion_status == "steigend":
        Speicherzustand[all_h_pos[-1]][2] = 0
        indices.reverse()
    n=0
    for i in indices:
        if inversion_status == "steigend":
            hPos_b = all_h_pos[i]
            hPos_r = all_h_pos[i+1]
        elif inversion_status == "fallend":
            hPos_b = all_h_pos[i+1]
            hPos_r = all_h_pos[i]

        Impuls = Speicherzustand[hPos_b][2]
        if Speicherzustand[hPos_b][3] > 0:
            f_imp_an_b = 0.1 * Speicherzustand[hPos_b][1] * speicher_param["A_Quer"]
        f_imp_an_b = min(f_imp_an_b, speicher_param["A_Quer"])
        counter = 0
        while Impuls > 0:
            theta_b = Speicherzustand[hPos_b][0]
            theta_r = Speicherzustand[hPos_r][0]


            V_b = Speicherzustand[hPos_b][1] * speicher_param["A_Quer"]
            dh_r = Speicherzustand[hPos_r][1]


            rho_b = __Modell_Stoffwerte("rho", theta_b)
            rho_r = __Modell_Stoffwerte("rho", theta_r)
            if inversion_status == "steigend":
                d_rho_inv = rho_b - rho_r
            elif inversion_status == "fallend":
                d_rho_inv = rho_r - rho_b
            
            d_V_an_b = f_imp_an_b * dh_r
            Speicherzustand[hPos_b][3] = Speicherzustand[hPos_b][3]\
                                         * V_b / (d_V_an_b + V_b)
            
            if Impuls**2 - g * d_rho_inv / rho_b * dh_r < 0:
                Impuls = 0
                Speicherzustand[hPos_b][2] = 0
                break
            Impuls_qdrt = (Impuls**2 * (1 
                                        - 2 * log((V_b * rho_b + d_V_an_b * rho_r)
                                                   / (V_b * rho_b)))
                                        - 2 * g * d_rho_inv / rho_b * dh_r)
            if Impuls_qdrt < 0:
                Impuls = 0
            else:
                Impuls = Impuls_qdrt**0.5
            
            theta_b = (V_b * theta_b * rho_b + d_V_an_b * theta_r * rho_r)\
                      / (V_b * rho_b + d_V_an_b * rho_r)
            V_b = (V_b * rho_b + d_V_an_b * rho_r)\
                  / __Modell_Stoffwerte("rho",theta_b)
            dh_r -= d_V_an_b / speicher_param["A_Quer"]

            Speicherzustand[hPos_r][1] = V_b / speicher_param["A_Quer"]
            Speicherzustand[hPos_b][1] = dh_r

            Speicherzustand[hPos_r][0] = theta_b
            Speicherzustand[hPos_b][0] = theta_r

            tmp = Speicherzustand[hPos_r][2]
            Speicherzustand[hPos_r][2] = Impuls
            Speicherzustand[hPos_b][2] = tmp

            tmp = Speicherzustand[hPos_r][3]
            Speicherzustand[hPos_r][3] = Speicherzustand[hPos_b][3]
            Speicherzustand[hPos_b][3] = tmp

            counter += 1

            if inversion_status == "steigend":
                if i + counter + 1 > len(Speicherzustand) - 1:
                    break
                hPos_b = all_h_pos[i+counter]
                hPos_r = all_h_pos[i+1+counter]
            elif inversion_status == "fallend":
                if i - counter < 0:
                    break
                hPos_b = all_h_pos[i+1-counter]
                hPos_r = all_h_pos[i-counter]
            n += 1
    Speicherzustand = __Modell_Aufraumen(Speicherzustand)
    return Speicherzustand

# // Funktion: Inversionsprüfung
# pruefen, ob Inversion vorliegt und ob fallende oder steigende dominiert
def __Modell_Inversionspruefung(Speicherzustand):
    # das wird immer wieder gemacht um sicherzustellen, dass die Reihenfolge
    # ueberall gleich bleibt

    all_h_pos = sorted(list(Speicherzustand))
    inversion = ""
    d_theta_Grenz = 1.0E-6
    inversionen = {
        "steigend" : 0,
        "fallend" : 0
    }
    for i in range(1, len(all_h_pos)-1):
        theta_i = Speicherzustand[all_h_pos[i]][0]
        theta_minus = Speicherzustand[all_h_pos[i-1]][0]
        theta_plus = Speicherzustand[all_h_pos[i+1]][0]
        inv_wert = theta_i - (theta_minus + theta_plus) / 2
        #inv_wert = round(inv_wert+0.0000001,4)
        if theta_i > theta_plus + d_theta_Grenz:
            if inv_wert > inversionen["steigend"]:
                inversionen["steigend"] = inv_wert
        elif theta_i < theta_minus - d_theta_Grenz:
            if inv_wert < inversionen["fallend"]:
                inversionen["fallend"] = inv_wert
    # Ende for
    # unteren Rand pruefen
    theta_plus = Speicherzustand[all_h_pos[1]][0]
    theta_i = Speicherzustand[all_h_pos[0]][0]
    inv_wert = theta_i - theta_plus
    if inv_wert > d_theta_Grenz:
        if inv_wert > inversionen["steigend"]:
            inversionen["steigend"] = inv_wert
    # oberen Rand pruefen
    theta_minus = Speicherzustand[all_h_pos[-2]][0]
    theta_i = Speicherzustand[all_h_pos[-1]][0]
    inv_wert = theta_i - theta_minus
    if inv_wert < -d_theta_Grenz:
        if inv_wert < inversionen["fallend"]:
            inversionen["fallend"] = inv_wert
    
    # nach groesstem Betrag schauen und damit festlegen ob steigende oder fallende
    # Inversion dominiert oder gar keine vorliegt
    if inversionen["fallend"] == 0 and inversionen["steigend"] == 0:
        inversion = "keine"
    elif abs(inversionen["steigend"]) > abs(inversionen["fallend"]):
        inversion = "steigend"
    elif abs(inversionen["steigend"]) <= abs(inversionen["fallend"]):
        inversion = "fallend"
    else:
        raise ValueError("Inversionen koennten nicht verglichen werden!")

    return inversion


# // Funktion: Inversionen auflösen
# Inversionen im Temperaturfeld aufloesen
def __Modell_Inversion(inversion_status, unten_oben, Vp_zu, dt, Speicherzustand):
    all_h_pos = sorted(list(Speicherzustand))
    any_inv = True
    all_inv_counter = 0

    while any_inv:

        any_inv = False
        indices = [i for i in range(0, (len(Speicherzustand)) - 1)]
        if inversion_status == "steigend":
            indices.reverse()
        for i in indices:
            V_b_lin_frac = 0
            hPos_r_next = -1
            d_theta_inv = 0
            d_rho_inv = 0

            inv = False
            BilanzMasse_1 = 0
            BilanzEnergie_1 = 0
            if inversion_status == "steigend":
                hPos_b = all_h_pos[i]
                hPos_r = all_h_pos[i+1]
                if i+2 <= len(all_h_pos)-1:
                    hPos_r_next = all_h_pos[i+2]
            elif inversion_status == "fallend":
                hPos_b = all_h_pos[i+1]
                hPos_r = all_h_pos[i]
                if i-1>=0:
                    hPos_r_next = all_h_pos[i-1]

            else:
                raise ValueError("'inversion_status' ist falsch") 
            theta_r = Speicherzustand[hPos_r][0]
            theta_b_kern = Speicherzustand[hPos_b][0]
            theta_b_grenz = theta_r
            if hPos_r_next > -1:
                theta_r_next = Speicherzustand[hPos_r_next][0]
            else:
                theta_r_next = None # B: ueberpruefen ob das hier eigentlich noetig ist
            if inversion_status == "steigend":
                d_theta_inv = theta_b_kern - theta_r
            else:
                d_theta_inv = theta_r - theta_b_kern
            if d_theta_inv > 1.0E-09: # ueberpruefen ob die zahl groesser sein darf
                any_inv = True
                inv = True
                if d_theta_inv > 1E-00:
                    BilanzMasse_1 = speicher_param["A_Quer"]\
                                    * sum([v[1] * __Modell_Stoffwerte("rho",v[0]) 
                                        for k,v in sorted(Speicherzustand.items())])
                    BilanzEnergie_1 = speicher_param["A_Quer"]\
                                     * sum([v[1] 
                                            * __Modell_Stoffwerte("rho",v[0]) 
                                            * __Modell_Stoffwerte("h",v[0]) 
                                        for k,v in sorted(Speicherzustand.items())])
            # Ende if d_theta_inv > 1e-09
            counter = 0
            while inv:
                inv = False
                """ Mischfaktor nimmt mit dem Volumenstrom zu, weil mit steigendem 
                Volumenstrom mehr Speicherquerschnitt vom aufsteigenden Fluid 
                eingenommen wird und daher mehr von den durchstroemten Zellen 
                eingemischt wird dadurch verlauft die Anpasung des bewegten Plugs 
                an die Umgebungstemperatur weniger abhaengig von Volumenstrom - 
                denn ein groesses bewegtes Plug braucht mehr Einmischung, um die 
                gleiche Temperaturentwickung entlang seiner Trajektorie zu nehmen
                $f_inv_an_b hat die Dimension m^3/m also Volumenzunahme der 
                Inversionsstroemung pro ueberwundener Speicherhoehe	
                """


                alt_d = 15.6
                alt_e = 0.165
                alt_f = 0.123
                alt_g = 1.08
                alt_h = 25
                alt_j = 5.54
                d = alt_d * 1
                e = alt_e * 1
                f = alt_f * 1
                g = alt_g * 1
                h = alt_h * 1
                j = alt_j * 1
                if unten_oben =="unten":
                    f_inv_an_b = (3 * Vp_zu - 0.04) * dt
                    f_inv_an_b = max(0.03 * dt, f_inv_an_b)
                else: 
                    
                    f_inv_an_b = (d * Vp_zu - e) * dt

                    f_inv_an_b = max(f * dt, f_inv_an_b)
                if Speicherzustand[hPos_b][3] == 0:
                    f_inv_an_b = f * dt
                if unten_oben == "unten":
                    f_inv_an_r = 3 - d_theta_inv / 10 - Vp_zu * 10
                else:
                    f_inv_an_r = g - d_theta_inv / h - Vp_zu * j
                f_inv_an_r = max(f_inv_an_r, 0)
                if Speicherzustand[hPos_b][3] == 0:
                    f_inv_an_r = 0.5



                V_b = Speicherzustand[hPos_b][1] * speicher_param["A_Quer"]
                # zelle aufloesen und raus aus dem <while inv>

                if V_b < 1.0E-9:
                    Speicherzustand[hPos_b][0] = Speicherzustand[hPos_r][0]
                    break

                V_b_lin = V_b * V_b_lin_frac
                V_b_kern = V_b - V_b_lin
                dh_r = Speicherzustand[hPos_r][1]
                
                if dh_r < 1.0E-12: # sehr kleine zellen aufloesen B: ACHTUNG ueberpruefen ob die zahl groesser sein koennte
                    dh_r_next = Speicherzustand[hPos_r_next][1] 
                    m_r = speicher_param["A_Quer"] *  dh_r * __Modell_Stoffwerte("rho",
                                                    Speicherzustand[hPos_r][0])
                    
                    m_r_next = speicher_param["A_Quer"] * dh_r_next * __Modell_Stoffwerte("rho", 
                                                    Speicherzustand[hPos_r_next][0])
                    E_r = m_r * __Modell_Stoffwerte("h", 
                                                    Speicherzustand[hPos_r][0])
                    E_r_next = m_r_next * __Modell_Stoffwerte("h", 
                                                    Speicherzustand[hPos_r_next][0])
                    theta_misch = __Modell_Stoffwerte("h_rev", ((E_r + E_r_next)
                                                                 / (m_r + m_r_next)))
                    rho_misch = __Modell_Stoffwerte("rho", theta_misch)

                    # B: ruhende Zelle mit der naechsten vermischen
                    Speicherzustand[hPos_r][1] = 0 # B: Zellenhoehe ruhende Zelle = 0
                    # B: naechste Zelle wird auch neu gebildet
                    Speicherzustand[hPos_r_next][0] = theta_misch
                    Speicherzustand[hPos_r_next][1] = (m_r + m_r_next) / rho_misch / speicher_param["A_Quer"] 
                    # B: Speicherzustand aufrauemen und die Funktion 'Inversion'
                    # verlassen. Verbleibende Inversionen werden spaeter aufgeloest
                    Speicherzustand = __Modell_Aufraumen(Speicherzustand)
                    return Speicherzustand
                # Ende if dh_r < 1.0E-12
                rho_b = __Modell_Stoffwerte("rho", Speicherzustand[hPos_b][0])
                rho_r = __Modell_Stoffwerte("rho", theta_r)
                if inversion_status == "steigend":
                    d_rho_inv = rho_b - rho_r
                elif inversion_status == "fallend":
                    d_rho_inv = rho_r - rho_b
                # B: volumen(hoehen)anteil der von ruhenden Zelle an die bewegte
                # abgegeben wird
                d_V_an_b = f_inv_an_b * dh_r
                # verhindern, dass mehr als das gesamte Zellvolumen aus der ruhen-
                # den Zelle entnommen wird
                d_V_an_b = min(d_V_an_b, dh_r * speicher_param["A_Quer"] * 0.99)
                # Impuls-Modell
                impulsQuad = Speicherzustand[hPos_b][2] ** 2 
                impulsQuad = (impulsQuad \
                              * (1 - 
                                 2 * log((V_b*rho_b + d_V_an_b*rho_r) / (V_b*rho_b))) 
                              - 2 * g * (d_rho_inv/rho_b) * dh_r)

                impulsQuad = max(0, impulsQuad)
                Speicherzustand[hPos_b][2] = impulsQuad**0.5

                Speicherzustand[hPos_b][3] = Speicherzustand[hPos_b][3] \
                                             * V_b / (d_V_an_b + V_b)
                A = (V_b_kern + V_b_lin / 2) * (theta_b_kern - theta_b_grenz)\
                    + d_V_an_b * (theta_r - theta_b_grenz)
                C = V_b_kern + V_b_lin + d_V_an_b # B: V_b_lin => linear mit der Temperatur zunehmndes Volumen d_V_an_b => durch Mischfaktor beeinfluste Volumenzunahme

                
                if  A / (theta_b_kern  - theta_b_grenz) >= C/2 : # Multiplikation ist rechnereisch guenstiger als Division?
                    V_b_kern = 2 * A / (theta_b_kern - theta_b_grenz) - C
                    V_b_lin = C - V_b_kern
                else:
                    V_b_kern = 0
                    V_b_lin = C
                    theta_b_kern = 2 * A / C + theta_b_grenz
                   
                # Zellhoehe (also Masse) bestimmen, die von der bewegten Zelle 
                # an die ruhende abgegeben wird
                d_V_an_r = 0
                if hPos_r_next > -1:
                    beta_rho = __Modell_Stoffwerte("beta_rho", theta_r) 
                    d_hPos_next = abs(hPos_r - hPos_r_next)
                    if inversion_status == "steigend":
                        theta_grenz_impuls = (theta_r_next 
                                              - f_inv_an_r 
                                                * impulsQuad 
                                                / (-2 * g * beta_rho * d_hPos_next ))
                    elif inversion_status == "fallend":
                        theta_grenz_impuls = (theta_r_next 
                                              - f_inv_an_r 
                                                * impulsQuad 
                                                / (2 * g * beta_rho * d_hPos_next ))
                    if theta_b_kern != theta_b_grenz:
                        d_V_an_r = (V_b_lin 
                                    * (theta_grenz_impuls - theta_b_grenz) 
                                    / (theta_b_kern - theta_b_grenz))
                    if ((theta_grenz_impuls > theta_b_kern 
                    or theta_grenz_impuls < theta_b_grenz)
                    and(inversion_status=="steigend")):
                        d_V_an_r = 0
                    if ((theta_grenz_impuls < theta_b_kern 
                    or theta_grenz_impuls > theta_b_grenz)
                    and(inversion_status=="fallend")):
                        d_V_an_r = 0
                # Ende if hPos_r_next > -1
                if d_V_an_r / speicher_param["A_Quer"] < 1E-06:
                    d_V_an_r = 0

                # neute Minimaltemperatur der bewegten  Zelle bestimmen
                d_theta_b_grenz = 0
                if d_V_an_r > 0:
                    d_theta_b_grenz = (theta_b_kern - theta_b_grenz)\
                                      * d_V_an_r / V_b_lin
                theta_b_grenz += d_theta_b_grenz

                # neue Zellparameter der ruhenden Zelle
                dh_r -= d_V_an_b / speicher_param["A_Quer"]
                theta_r = ((theta_r * dh_r * speicher_param["A_Quer"]
                            + (theta_b_grenz - 0.5 * d_theta_b_grenz) * d_V_an_r)
                           / (dh_r * speicher_param["A_Quer"] + d_V_an_r))
                dh_r += d_V_an_r / speicher_param["A_Quer"]
                V_b_lin -= d_V_an_r

                V_b = V_b_lin + V_b_kern 
                if V_b == 0: # ACHTUNG  muss es "==" sein? (<= )
                    break
                V_b_lin_frac = V_b_lin / V_b
                
                theta_b_misch = ((V_b_kern * theta_b_kern 
                                  + (theta_b_kern + theta_b_grenz) / 2 * V_b_lin)
                                 / V_b)
                
                

                Speicherzustand[hPos_r][1] = V_b / speicher_param["A_Quer"]
                Speicherzustand[hPos_b][1] = dh_r
                Speicherzustand[hPos_r][0] = theta_b_misch
                Speicherzustand[hPos_b][0] = theta_r
                tmp = Speicherzustand[hPos_r][2]
                Speicherzustand[hPos_r][2] = Speicherzustand[hPos_b][2]
                Speicherzustand[hPos_b][2] = tmp
                tmp = Speicherzustand[hPos_r][3]
                Speicherzustand[hPos_r][3] = Speicherzustand[hPos_b][3]
                Speicherzustand[hPos_b][3] = tmp

                # Ueberpruefung des naechsten Zellpaares vorbereiten
                counter += 1
                hPos_r_next = -1
                if inversion_status == "steigend":
                    if i + counter + 2 > len(Speicherzustand) - 1: 
                        break
                    hPos_b = all_h_pos[i + counter]
                    hPos_r = all_h_pos[i + counter + 1]
                    if i + counter + 2 <= len(all_h_pos) - 1: 
                        hPos_r_next = all_h_pos[i + counter + 2]
                # Ende if invesion_status=="steigend"
                if inversion_status == "fallend":
                    if i - counter - 1 < 0:
                        break
                    hPos_b = all_h_pos[i - counter + 1]
                    hPos_r = all_h_pos[i - counter]
                    if i - counter - 1 >= 0: 
                        hPos_r_next = all_h_pos[i - counter - 1]
                # Ende elif inversion_status == "fallend"

                theta_r = Speicherzustand[hPos_r][0]
                if hPos_r_next > -1: 
                    theta_r_next = Speicherzustand[hPos_r_next][0]
                if inversion_status == "steigend":
                    d_theta_inv = Speicherzustand[hPos_b][0] - theta_r
                if inversion_status == "fallend":
                    d_theta_inv = theta_r - Speicherzustand[hPos_b][0]
                inv = True if d_theta_inv > 0 else False

            # Ende <while inv>
            if BilanzMasse_1 > 0:
                BilanzMasse_2  = speicher_param["A_Quer"]\
                                 * sum([v[1] * __Modell_Stoffwerte("rho",v[0]) 
                                    for k,v in sorted(Speicherzustand.items())])
                BilanzEnergie_2 = speicher_param["A_Quer"]\
                                * sum([v[1] 
                                        * __Modell_Stoffwerte("rho",v[0]) 
                                        * __Modell_Stoffwerte("h",v[0]) 
                                        for k,v in sorted(Speicherzustand.items())])
                m_alt = speicher_param["A_Quer"] * Speicherzustand[hPos_b][1]\
                        * __Modell_Stoffwerte("rho", Speicherzustand[hPos_b][0])
                E_alt = m_alt * __Modell_Stoffwerte("h", Speicherzustand[hPos_b][0])
                m_neu = m_alt - (BilanzMasse_2 - BilanzMasse_1)
                E_neu = E_alt - (BilanzEnergie_2 - BilanzEnergie_1)
                theta_neu = __Modell_Stoffwerte("h_rev", E_neu / m_neu)
                rho_neu = __Modell_Stoffwerte("rho", theta_neu)
                Speicherzustand[hPos_b][0] = theta_neu
                Speicherzustand[hPos_b][1] = m_neu / speicher_param["A_Quer"]\
                                             / rho_neu
            # Ende if Bilanzmasse_1 > 1
        # Ende for i in indices
        all_inv_counter += 1
    # Ende <while any_inv>
    Speicherzustand = __Modell_Aufraumen(Speicherzustand)
    

    return Speicherzustand
        
# // Funktion: Zustrom
def __Modell_Zustrom(unten_oben, aktuellSekunden, ausgabezeit, Begleitdaten,
                     Zeitabstand, Ausgabewerte, Speicherzustand):
    """
    Erzeugung neuer Wasserzellen, also das eintrende Wasser 
    Es werden die Position und Eigenschaften der Zell berechnet.
    
    """
    all_h_pos = sorted(list(Speicherzustand))
    h_WS = max(all_h_pos) + Speicherzustand[max(all_h_pos)][1] / 2
    if unten_oben == "unten":
        h_zu_min = speicher_param["H_B_UK_Dif"] + 0.01 * speicher_param["H_RS_Dif"]
        h_zu_max = speicher_param["H_B_UK_Dif"] + 0.99 * speicher_param["H_RS_Dif"]
        theta_VL_RL = "T_RL"
        Ausgabename_h_zu = "h_zu_unten"
    elif unten_oben == "oben":
        # Hoehenposition der Wassersaeule
        h_zu_min = h_WS - (speicher_param["H_WS_OK_Dif"] 
                           + 0.99 * speicher_param["H_RS_Dif"])
        h_zu_max = h_WS - (speicher_param["H_WS_OK_Dif"] 
                           + 0.01 * speicher_param["H_RS_Dif"])
        theta_VL_RL = "T_VL"

    else:
        raise ValueError("'unten' oder 'oben' spezifizieren")
    Ausgabewerte[theta_VL_RL] = {}

    h_zu_theta = (h_zu_max + h_zu_min) / 2
    d_theta_h_zu_theta = 1.0E12

    theta_zu = Begleitdaten[aktuellSekunden][theta_VL_RL]
    m_punkt_zu = Begleitdaten[aktuellSekunden]["m_Punkt"]
    rho_zu = __Modell_Stoffwerte("rho", theta_zu)
    Vp_zu = m_punkt_zu / rho_zu  

    dh_zu = Vp_zu * Zeitabstand / speicher_param["A_Quer"] #meter
    v_zu = Vp_zu / (pi * 2 * speicher_param["R_Dif"] * speicher_param["H_RS_Dif"])  # Austrittsgeschwindigkeit bei perfekter Spaltfuellung in m/s

    # das Plug im Hoehenbereich bestimmen, dessen Temperatur am naechsten an der Zu-
    # trittstemperatur liegt, diesem dann das Fluid zufuehren
    if unten_oben == "unten":
        for hPos in all_h_pos:
            if hPos < h_zu_min:
                continue
            elif (hPos <= h_zu_max and hPos >= h_zu_min):
                if abs(theta_zu - Speicherzustand[hPos][0]) < d_theta_h_zu_theta:
                    h_zu_theta = hPos
                    d_theta_h_zu_theta = abs(theta_zu - Speicherzustand[hPos][0])
            elif hPos > h_zu_max:
                break
            else:
                raise ValueError("hPos comparison error in __Modell_Zustrom()")
    # Ende if unten_oben == "unten"
    elif unten_oben == "oben":
        for hPos in reversed(all_h_pos):
            if hPos > h_zu_max:
                continue
            elif (hPos <= h_zu_max and hPos >= h_zu_min):
                if abs(theta_zu - Speicherzustand[hPos][0]) < d_theta_h_zu_theta:
                    h_zu_theta = hPos
                    d_theta_h_zu_theta = abs(theta_zu - Speicherzustand[hPos][0])
            elif hPos < h_zu_min:
                break
            else:
                raise ValueError("hPos Comparison error in __Modell_Zustrom()")
    # Ende elif unten_oben == "oben"

    fak_neben = 0


    # extra Genauigkeit mit den richtigen Angaben, die aber sehr schwer zu bekommen sind, also fragwürdig ob sinnvoll
    if speicher_param["nebenstrom"] == 1:
        if unten_oben == "oben":
            dh_zu, theta_zu, fak_neben, Speicherzustand = __Modell_Nebenstrom_zu(theta_zu, Vp_zu, h_WS, dh_zu, Speicherzustand)
    else:
        fak_neben = 0
    index_Nachbar = __Modell_find_index_h_pos(h_zu_theta, all_h_pos)
    hPos_Nachbar = all_h_pos[index_Nachbar]
    #print(f"hpos nachbar:{hPos_Nachbar}")

    if theta_zu >= Speicherzustand[hPos_Nachbar][0]:
        fallend_steigend = "steigend"

    elif theta_zu < Speicherzustand[hPos_Nachbar][0]:
        fallend_steigend = "fallend"

    else:
        raise ValueError("fallend_steigend konnte nicht erzeugt werde")
    #print(f"H zu zustrom before: {h_zu}")
    f_v_max = 1
    v_zu_effektiv = v_zu * f_v_max * (1+fak_neben)
    #hPos_r = all_h_pos[index_Nachbar]
    Impuls = v_zu
    counter = 0
    # unklar ob das verbessert oder verschlechtert, führt manchmal zu komischen Sachen
    # deshalb auskommentiert (sehr präzise und genaue h_zu ist prinzipiell nicht so wichtig, wir reden da um cm)
    # keine große Auswirkung auf das Temperaturprofil, wenn man 40m Speicher betrachtet. 
    # if fallend_steigend == "steigend" and unten_oben == "unten":
    #     while Impuls > 0:
    #         dh_Nachbar = Speicherzustand[hPos_Nachbar][1]
    #         d_rho_inv = rho_zu - __Modell_Stoffwerte("rho",
    #                                              Speicherzustand[hPos_Nachbar][0])
    #         Impuls_quadrat = (0.5 * 2 * g * d_rho_inv / rho_zu * dh_Nachbar 
    #                           + Impuls**2)
    #         if Impuls_quadrat < 0:
    #             break
    #         Impuls = Impuls_quadrat**0.5
    #         counter += 1
    #         if index_Nachbar - counter <= 0:
    #             break
    #         hPos_Nachbar = all_h_pos[index_Nachbar - counter]
    #     # Ende while
    #     h_zu = hPos_Nachbar - Speicherzustand[hPos_Nachbar][1] / 2
    # # Ende if
    # elif fallend_steigend == "fallend" and unten_oben == "oben":
    #     while Impuls > 0:
    #         dh_Nachbar = Speicherzustand[hPos_Nachbar][1]
    #         d_rho_inv = - rho_zu + __Modell_Stoffwerte("rho",
    #                                              Speicherzustand[hPos_Nachbar][0])
    #         Impuls_quadrat = 2 * g * d_rho_inv / rho_zu * dh_Nachbar + Impuls**2
    #         if Impuls_quadrat < 0:
    #             break
    #         Impuls = Impuls_quadrat**0.5
    #         counter += 1
    #         if (index_Nachbar + counter) >= (len(all_h_pos) - 1):     # len(lst) - 1 => letztes index der Liste
    #             break
    #         hPos_Nachbar = all_h_pos[index_Nachbar + counter]
    #     # Ende while
    #     h_zu = hPos_Nachbar + Speicherzustand[hPos_Nachbar][1] / 2
    # # # Ende elif
    # else:

    h_zu = h_zu_theta
    if theta_zu >= Speicherzustand[h_zu_theta][0]:
        h_zu += Speicherzustand[h_zu_theta][1] / 2
    elif theta_zu < Speicherzustand[h_zu_theta][0]:
        h_zu -= Speicherzustand[h_zu_theta][1] / 2
    else:
        raise ValueError("theta_zu koennte nicht verglichen werden")

    # Ende else
    if dh_zu > 0:
        Speicherzustand[h_zu] = [theta_zu, dh_zu, 0, v_zu_effektiv]
    Ausgabewerte[theta_VL_RL][ausgabezeit] = theta_zu

    return Speicherzustand

# // Funktion: Abstrom
def __Modell_Abstrom(unten_oben, aktuellSekunden, ausgabezeit, Ausgabewerte,
                     Begleitdaten, Zeitabstand, Speicherzustand):
    """
    Entfernen der Wasserzellen ist hier. Es werden die Zellen im Bereich des Diffusors
    zusammengemischt und entfernt. Es wird die Enthalpie der gemischten Zellen berechnet
    und anschließend die Ausstromtemperatur berechnet
    """
    Speicherzustand = __Modell_Aufraumen(Speicherzustand)
    all_h_pos = sorted(list(Speicherzustand))
    lastKey = max(all_h_pos)
    h_WS = lastKey + Speicherzustand[lastKey][1] / 2
    theta_diffusor = __Modell_Temperatur_Diffusorhoehe(unten_oben, h_WS, Speicherzustand)
    F_ab = Begleitdaten[aktuellSekunden]["m_Punkt"] / __Modell_Stoffwerte("rho", theta_diffusor)
    if unten_oben == "unten":
        h_ab_min = speicher_param["H_B_UK_Dif"]
        h_ab_max = speicher_param["H_B_UK_Dif"] + speicher_param["H_RS_Dif"]
        Ausgabename = "T_RL"
    elif unten_oben == "oben":
        h_ab_min = h_WS - speicher_param["H_WS_OK_Dif"] - speicher_param["H_RS_Dif"]
        h_ab_max = h_WS - speicher_param["H_WS_OK_Dif"] # Änderung um die Abweichungen zu beheben (oben beim Entladen wenn der Speicher fast leer ist)
        Ausgabename = "T_VL"

    F_neben = 0
    theta_neben = 0

    # TODO FIXME
    if unten_oben == "oben":
        theta_neben, F_neben, Speicherzustand = __Modell_Nebenstrom_ab(F_ab,h_WS, Zeitabstand,Speicherzustand)

    Ausgabewerte[Ausgabename] = {} 
    dh_ab_rel= {}
    rho_mittel = 0

    if unten_oben == "unten":
        for hPos in all_h_pos:
            theta = Speicherzustand[hPos][0]
            h_zelle_min = hPos - Speicherzustand[hPos][1] / 2
            h_zelle_max = hPos + Speicherzustand[hPos][1] / 2
            if h_zelle_max < h_ab_min:
                continue

            elif h_zelle_min <= h_ab_min and h_zelle_max >= h_ab_max:
                dh_ab_rel[hPos] = 1
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho",theta)

            elif h_zelle_min <= h_ab_min and h_zelle_max >= h_ab_min:
                dh_ab_rel[hPos] = (h_zelle_max - h_ab_min) / speicher_param["H_RS_Dif"]
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho", theta)

            elif h_zelle_min >= h_ab_min and h_zelle_max <= h_ab_max:
                dh_ab_rel[hPos] = Speicherzustand[hPos][1] / speicher_param["H_RS_Dif"]
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho", theta)

            elif h_zelle_min <= h_ab_max and h_zelle_max >= h_ab_max:
                dh_ab_rel[hPos] = (h_ab_max - h_zelle_min) / speicher_param["H_RS_Dif"]
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho", theta)

            elif h_zelle_min > h_ab_max:
                break
        # ende for hpos
    # ende if
    elif unten_oben == "oben":
        for hPos in reversed(all_h_pos):
            theta = Speicherzustand[hPos][0]
            h_zelle_min = hPos - Speicherzustand[hPos][1] / 2
            h_zelle_max = hPos + Speicherzustand[hPos][1] / 2
            if h_zelle_min > h_ab_max:
                continue

            elif h_zelle_min <= h_ab_min and h_zelle_max >= h_ab_max:
                dh_ab_rel[hPos] = 1
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho",theta)

            elif h_zelle_min <= h_ab_min and h_zelle_max >= h_ab_min:
                dh_ab_rel[hPos] = (h_zelle_max - h_ab_min) / (speicher_param["H_RS_Dif"])
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho", theta)

            elif h_zelle_min >= h_ab_min and h_zelle_max <= h_ab_max:
                dh_ab_rel[hPos] = Speicherzustand[hPos][1] / (speicher_param["H_RS_Dif"])
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho", theta)

            elif h_zelle_min <= h_ab_max and h_zelle_max >= h_ab_max:
                dh_ab_rel[hPos] = (h_ab_max - h_zelle_min) / (speicher_param["H_RS_Dif"])
                rho_mittel += dh_ab_rel[hPos] * __Modell_Stoffwerte("rho", theta)

            elif h_zelle_max < h_ab_min:
                break
        # ende for hpos
    # ende elif
    sum_dh_ab_rel = sum([v for v in dh_ab_rel.values()])
    if abs(sum_dh_ab_rel - 1) > 1.0E-06:
        raise ValueError("sum_dh_ab_rel ist nicht gleich 1 : %s" %sum_dh_ab_rel)

    m_ab = 0
    H_ab = 0
    theta_ab = 0

    if F_neben > 0:
        rho_neben = __Modell_Stoffwerte("rho", theta_neben)
        h_neben = __Modell_Stoffwerte("h", theta_neben)
        m_ab = rho_neben * F_neben * Zeitabstand
        
        H_ab = m_ab * h_neben

    V_ab_soll = Begleitdaten[aktuellSekunden]["m_Punkt"] / rho_mittel * Zeitabstand
    V_ab_ist = F_neben * Zeitabstand /3600

    while (abs(V_ab_soll - V_ab_ist) > 1.0E-06):
        dh_ab_ges = (V_ab_soll  - V_ab_ist) / speicher_param["A_Quer"]
        for hPos in sorted(list(dh_ab_rel)):
            theta = Speicherzustand[hPos][0]
            dh_ab = dh_ab_rel[hPos] * dh_ab_ges
            if dh_ab > Speicherzustand[hPos][1]:
                raise ValueError("Aus der Zelle auf pos. %s soll zu viel entnommen werden. dh_ab = %s, rel = %s, dh_ab_ges = %s"%(hPos, dh_ab, dh_ab_rel[hPos], dh_ab_ges))
            Speicherzustand[hPos][1] -= dh_ab

            H_ab  += (dh_ab * speicher_param["A_Quer"] 
                      * __Modell_Stoffwerte("rho", theta)
                      * __Modell_Stoffwerte("h", theta)
                      )
            m_ab += dh_ab * speicher_param["A_Quer"]\
                    * __Modell_Stoffwerte("rho",theta)
        # ende for hpos

        theta_ab = __Modell_Stoffwerte("h_rev", H_ab / m_ab)
        rho_ab = __Modell_Stoffwerte("rho", theta_ab)
        V_ab_ist = m_ab/rho_ab
    # ende while
    Ausgabewerte[Ausgabename][ausgabezeit] = theta_ab
    Speicherzustand = __Modell_Aufraumen(Speicherzustand)

    return Speicherzustand


# // Funktion: Wärmeleitung
def __Modell_Waermeleitung(Zeitabstand, thetaRand, q_punkt, 
                           Fundamentzustand, Speicherzustand):
    """
    Berechnet die Wärmeleitung im Inneren des Speichers anhand der Wärmeleitungsgleichung

    """

    Gesamtzustand = {}
    Gesamtzustand = {**Fundamentzustand,**Speicherzustand}          # Fundament- und Speicherzustand werden in eine Variable geschrieben
    all_h_pos = (sorted(list(Gesamtzustand)))                       # Sortierung aller Höhenpositionen in einer Liste
    all_h_pos.reverse()                                             # Reihenfolge vertauschen, Index 0 ist oben, Speicher wird von oben nach unten durchgegangen
    thetaWL = [Gesamtzustand[h][0] for h in all_h_pos]              # Liste aller Temperaturen der Zellen
    dx = [Gesamtzustand[h][1] for h in all_h_pos]                   # Liste aller Dicken der Zellen
    maxIndex = len(thetaWL) - 1                                     # Index des letzten Eintrags der Liste
    dx_ = [None for h in all_h_pos]                                 # Erstellen von 7 Listen mit der Länge all_h_pos und dem Inhalt "None"
    tlf_mod = [None for h in all_h_pos]
    a_WL = [None for h in all_h_pos]
    b_WL = [None for h in all_h_pos]
    c_WL = [None for h in all_h_pos]
    d_WL = [None for h in all_h_pos]
    lambda_ = [None for h in all_h_pos]

    for j in range(0, maxIndex):
        dx_[j] = (dx[j] + dx[j+1]) / 2                              # arithmetisches Mittel der Dicke von der Zelle und der Zelle darunter (weil von oben nach unten gegangen wird) bis zum vorletzten Element der Liste
    for j in range(0, len(thetaWL)):
    #for j in range(0, maxIndex):
        if all_h_pos[j] > 0:                                        # Dichte und cp für Zellen im Speichermedium und im Fundament berechnen
            rho = __Modell_Stoffwerte("rho", thetaWL[j])            # für Dichte die Temperatur der Zelle selbst nehmen
            cp = __Modell_Stoffwerte("cp", thetaWL[j+1])            # für Wärmekapazität die Temperatur der Zelle darunter nehmen
        else:
            rho = __Modell_Stoffwerte("rho_Fundament")
            cp = __Modell_Stoffwerte("cp_Fundament")
        tlf_mod[j] = Zeitabstand / (rho * cp * 2 * dx[j])          
    # ende for
    
    for j in range(0, len(thetaWL)-1):                                          # Wärmeleitfähigkeit für alle Zellen im Speicher und Fundament berechnen
        if all_h_pos[j] > 0:
            Lambda = __Modell_Stoffwerte("lambda",thetaWL[j])
        if all_h_pos[j+1] > 0:
            Lambda_plus = __Modell_Stoffwerte("lambda", thetaWL[j+1])
        if all_h_pos[j] < 0:
            Lambda = __Modell_Stoffwerte("lambda_Fundament")
        if all_h_pos[j+1] < 0:
            Lambda_plus = __Modell_Stoffwerte("lambda_Fundament")
        lambda_[j] = (dx[j] + dx[j+1]) / (dx[j]/Lambda + dx[j+1]/Lambda_plus)   # mittlere Wärmeleitfähigkeit?? der Zelle und der Zelle darunter
    # ende for

    lambda_0 = __Modell_Stoffwerte("lambda", thetaWL[0])
    c_WL[0] = -tlf_mod[0] * lambda_[0] / dx_[0]
    b_WL[0] = 1 + tlf_mod[0] * (lambda_[0] / dx_[0] + lambda_0 / (dx[0] / 2))
    d_WL[0] = thetaWL[0] * (1 - tlf_mod[0] * (lambda_[0] / dx_[0] + lambda_0 / (dx[0] / 2)))\
                + thetaWL[1] * tlf_mod[0] * lambda_[0] / dx_[0]\
                + 2 * thetaRand * tlf_mod[0] * lambda_0 / (dx[0] / 2)
    rho_f = __Modell_Stoffwerte("rho_Fundament")
    cp_f = __Modell_Stoffwerte("cp_Fundament")
    a_WL[maxIndex] = -tlf_mod[maxIndex] * lambda_[maxIndex - 1] / dx_[maxIndex - 1]
    b_WL[maxIndex] = 1 + tlf_mod[maxIndex] * lambda_[maxIndex - 1] / dx_[maxIndex - 1]
    d_WL[maxIndex] = thetaWL[maxIndex] * (1 - tlf_mod[maxIndex] * lambda_[maxIndex - 1] / dx_[maxIndex - 1])\
                    + thetaWL[maxIndex - 1] * tlf_mod[maxIndex] * lambda_[maxIndex - 1] / dx_[maxIndex - 1]\
                    - q_punkt * Zeitabstand / (dx[maxIndex] * rho_f * cp_f)


    for j in range(1,maxIndex):                                                 # Berechnung der Koeffizienten zur Lösung des linearen Gleichungssystems der Wärmeleitungsgleichung für jede Zelle
        a_WL[j] = -tlf_mod[j] * lambda_[j-1] / dx_[j-1]
        b_WL[j] = 1 + tlf_mod[j] * (lambda_[j] / dx_[j] + lambda_[j-1] / dx_[j-1])
        c_WL[j] = -tlf_mod[j] * lambda_[j] / dx_[j]
        d_WL[j] = thetaWL[j] * (1 - tlf_mod[j] * (lambda_[j-1] / dx_[j-1] + lambda_[j] / dx_[j]))\
                 + thetaWL[j+1] * tlf_mod[j] * lambda_[j] / dx_[j]\
                 + thetaWL[j-1] * tlf_mod[j] * lambda_[j-1] / dx_[j-1]

    thetaWL = __Modell_TDMASolve(a_WL, b_WL, c_WL, d_WL)                        # Berechnung der neuen Temperatur jeder Zelle durch Lösung der Wärmeleitungsgleichung
    j = 0
    thetaWL[0] = thetaRand
    for hPos in all_h_pos:                                                      # neue Temperaturen und Zellhöhen in Speicherzustand und Fundamentzustand schreiben
        if hPos < 0: 
            Fundamentzustand[hPos][0] = thetaWL[j]
        if hPos > 0:
            masse = Speicherzustand[hPos][1]\
                    * __Modell_Stoffwerte("rho", Speicherzustand[hPos][0])
            Speicherzustand[hPos][1] = masse / __Modell_Stoffwerte("rho", thetaWL[j])
            Speicherzustand[hPos][0] = thetaWL[j]
        j += 1                                                                  # Frage: hPos wird nicht korrigiert, wann passiert das? in Modell_Aufräumen!
            
    return Fundamentzustand, Speicherzustand

# // Funktion: TDMASolve
def __Modell_TDMASolve(aL,bL, cL, dL):                              # Tridiagonales Matrixalgorithmus-Verfahren
    """
    Löst das lineare Gleichungssystem der Wärmeleitungsgleichung

    """
    x = [None for i in dL]
    n = len(dL) - 1 # max index
    for i in range(1, n + 1): # n muss inclusive sein
        m = aL[i] / bL[i-1]
        bL[i] -= m * cL[i-1]
        dL[i] -= m * dL[i-1]
    x[n] = dL[n] / bL[n]

    for i in reversed(range(0, n)):
        x[i] = (dL[i] - cL[i] * x[i+1]) / bL[i]                     # Berechnung der neuen Temperatur nach Durchführung der Wärmeleitung

    return x

# // Funktion: Kapazitäten
def __Modell_Kapazitaeten(dt_sub, T_amb, Kapazitaeten, Speicherzustand):
    """
    Berechnet den Wärmeübergang zwischen Mantel und Wasser und anschließend zwischen Mantel und Umgebung
    """
    U_Wert = 800 # Wärmeübergangskoeffizient !!!! (eigentlich) ALPHA  (von wasser zu innere Manteloberfläche) W/m2K (nicht wie in config)
    E_an_W = {} # von bauteilkapa an wasser energie
    all_h_pos_K = sorted(list(Kapazitaeten))
    all_h_pos_W = sorted(list(Speicherzustand))
    Fuellstand = all_h_pos_W[-1] + Speicherzustand[all_h_pos_W[-1]][1] / 2          # Position der obersten Zelle + halbe Zellhöhe


    if Fuellstand > speicher_param["H_Mantel"]:
        raise ValueError("Wasserspiegel hoeher als Speicher")

    counter_Kap = 0
    for hPosW in all_h_pos_W:
        dh_W = Speicherzustand[hPosW][1]                                # Position der Zellen im Wasser (= Wasserzelle)
        uG_W = hPosW - dh_W / 2                                         # Unterkante der Zelle
        oG_W = hPosW + dh_W / 2                                         # Oberkante der Zelle
        theta_W = Speicherzustand[hPosW][0]                             # Temperatur der Zelle
        E_an_W[hPosW] = 0
        next_K = True
        while next_K:
            next_K = False
            hPosK = all_h_pos_K[counter_Kap]                            # Position der Zelle im Mantel (= Mantelzelle)
            dh_K = Kapazitaeten[hPosK][1]                               # Höhe der Mantelzelle, die in Kapazitäten definiert ist
            uG_K = hPosK - dh_K / 2                                     # Unterkante dieser Zelle
            oG_K = hPosK + dh_K / 2                                     # Oberkante dieser Zelle
            theta_K = Kapazitaeten[hPosK][0]                            # Temperatur der Mantelzelle
            A = 0
            if uG_W >= uG_K and oG_W < oG_K:                            # Wasserzelle liegt in der Mantelzelle, die gerade betrachtet wird
                A = dh_W * pi * 2 * speicher_param["R_innen"]               # es wird Höhe der Wasserzelle zur Berechnung der Kontaktfläche der Zelle zwischen Wasser und Mantel verwendet
            elif uG_W >= uG_K and oG_W >= oG_K:                         # Wasserzelle liegt unten in der Mantelzelle, geht aber oben über diese hinaus
                A = (oG_K - uG_W) * pi * 2 * speicher_param["R_innen"]      # es wird nur der Bereich in dieser Mantelzelle zur Berechnung verwendet und dann erneut in die Schleife gegangen mit der nächsten Mantelzelle
                next_K = True
                counter_Kap += 1
            elif uG_W < uG_K and oG_W < oG_K:                           # Wasserzelle liegt unterhalb der Mantelzelle, aber geht oben nicht über diese hinaus
                A = (oG_W - uG_K) * pi * 2 * speicher_param["R_innen"]      # es wird nur der Bereich in dieser Mantelzelle zur Berechnung verwendet
            elif uG_W < uG_K and oG_W >= oG_K:                          # Wasserzelle liegt unterhalb der Mantelzelle und geht bis zur Oberkante der Zelle oder darüber hinaus
                A = dh_K * pi * 2  * speicher_param["R_innen"]              # es wird die Höhe der Mantelzelle verwendet und dann erneut in die Schleife gegangen mit der nächsten Mantelzelle
                next_K = True
                counter_Kap += 1


            E = dt_sub * speicher_param["alpha_water_innerwall"] * A * (theta_K - theta_W)               # Berechnung wie viel Energie zwischen Mantel und Wasser übertragen wird
            E_an_W[hPosW] += E                                          # übertragene Energie an jede Wasserzelle wird berechnet
            Kapazitaeten[hPosK][0] -= E / Kapazitaeten[hPosK][2]        # Temperatur des Mantels wird korrigiert mit übertragener Energie
        # ende while next_k

    # ende hposw
    for hPosW in all_h_pos_W:
        theta_W = Speicherzustand[hPosW][0]                                             # Temperatur der Wasserzelle
        dh_W = Speicherzustand[hPosW][1]                                                # Höhe der Wasserzelle
        m_W = __Modell_Stoffwerte("rho", theta_W) * dh_W * speicher_param["A_Quer"]     # Masse der Wasserzelle
        H = m_W * __Modell_Stoffwerte("h", theta_W) + E_an_W[hPosW]                     # neue Enthalpie der Wasserzelle
        theta_W_neu = __Modell_Stoffwerte("h_rev", H/m_W)                               # neue Temperatur der Wasserzelle
        rho_W_neu = __Modell_Stoffwerte("rho", theta_W_neu)                             # neue Dichte der Wasserzelle
        Speicherzustand[hPosW][0] = theta_W_neu                                         # neue Temperatur in Speicherzustand schreiben
        Speicherzustand[hPosW][1] = m_W / rho_W_neu / speicher_param["A_Quer"]          # neue Höhe der Zelle berechnen und in Speicherzustand schreiben
    E_Verlust_Mantel = 0
    for hPosK in all_h_pos_K:
        if hPosK + Kapazitaeten[hPosK][1]/2 > Fuellstand:                               # nur für Zellen unter dem Wasserspiegel durchführen
            continue
        E_Verlust = (dt_sub * (Kapazitaeten[hPosK][0] - T_amb)
                     * speicher_param["U_Mantel"]
                     * pi * 2 * speicher_param["R_innen"]
                     * Kapazitaeten[hPosK][1])                                          # Verlust des Speichermantels an die Umgebung berechnen, Frage: Warum R_innen und nicht der Radius außen am Speicher?
        Kapazitaeten[hPosK][0] -= E_Verlust / Kapazitaeten[hPosK][2]                    # neue Temperatur der Mantelzelle berechnen
        E_Verlust_Mantel += E_Verlust                                                   # Gesamtverlust der Energie über den Mantel an die Umgebung
    return E_Verlust_Mantel, Kapazitaeten, Speicherzustand                              # Frage: Müsste man nicht Kapazitäten[2] für die neue Temperatur berechnen, da sich cp und rho ändern?

# // Funktion: Temperatur Diffusorhöhe
def __Modell_Temperatur_Diffusorhoehe(unten_oben, h_WS, Speicherzustand):

    dh_rel = {}
    theta_mittel = 0
    all_h_pos = sorted(list(Speicherzustand))
    if unten_oben == "unten":
        h_ab_min = speicher_param["H_B_UK_Dif"]
        h_ab_max = speicher_param["H_B_UK_Dif"] + speicher_param["H_RS_Dif"]
        for hPos in all_h_pos:
            h_plug_min = hPos - Speicherzustand[hPos][1] / 2
            h_plug_max = hPos + Speicherzustand[hPos][1] / 2

            if h_plug_max < h_ab_min:
                continue
            elif (h_plug_min <= h_ab_min and h_plug_max >= h_ab_max):
                dh_rel[hPos] = 1
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]
            elif (h_plug_min <= h_ab_min and h_plug_max >= h_ab_min):
                dh_rel[hPos] = (h_plug_max - h_ab_min) / speicher_param["H_RS_Dif"]
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]
            elif (h_plug_min >= h_ab_min and h_plug_max <= h_ab_max):
                dh_rel[hPos] = Speicherzustand[hPos][1] / speicher_param["H_RS_Dif"]
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]
            elif (h_plug_min <= h_ab_max and h_plug_max >= h_ab_max):
                dh_rel[hPos] = (h_ab_max - h_plug_min) / speicher_param["H_RS_Dif"]
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]
            elif (h_plug_min > h_ab_max):
                break
    elif unten_oben == "oben":
        h_ab_min = h_WS - speicher_param["H_WS_OK_Dif"] - speicher_param["H_RS_Dif"]
        h_ab_max = h_WS - speicher_param["H_WS_OK_Dif"]
        for hPos in reversed(all_h_pos):
            h_plug_min = hPos - Speicherzustand[hPos][1] / 2
            h_plug_max = hPos + Speicherzustand[hPos][1] / 2

            if h_plug_min > h_ab_max:
                continue

            elif (h_plug_min <= h_ab_min and h_plug_max >= h_ab_max):
                dh_rel[hPos] = 1
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]

            elif (h_plug_max >= h_ab_max and h_plug_min <= h_ab_max):
                dh_rel[hPos] = (h_ab_max - h_plug_min) / speicher_param["H_RS_Dif"]
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]

            elif (h_plug_min >= h_ab_min and h_plug_max <= h_ab_max):
                dh_rel[hPos] = Speicherzustand[hPos][1] / speicher_param["H_RS_Dif"]
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]

            elif (h_plug_min <= h_ab_min and h_plug_max >= h_ab_min):
                dh_rel[hPos] = (h_plug_max - h_ab_min) / speicher_param["H_RS_Dif"]
                theta_mittel += dh_rel[hPos] * Speicherzustand[hPos][0]

            elif h_plug_max < h_ab_min:
                break
    # if abs(sum(dh_rel.values()) - 1)>1.0E-06:
    #     print("ERROR IN TEMPERATUR DIFFUSORHOEHOE %s" %(sum(dh_rel.values())))
        #raise ValueError("dh_rel muss 1 ergeben (Diffusor temperatur) %s"%(sum(dh_rel)))
    return theta_mittel
    
# // Funktion: Nebenstrom zu
# \\ nur für Speicher Dessau notwendig
def __Modell_Nebenstrom_zu(theta_zu, F_zu, h_WS, dh_zu, Speicherzustand):

    all_h_pos = sorted(list(Speicherzustand))
    fak_neben = 0
    rho_1 = __Modell_Stoffwerte("rho", theta_zu)
    u_1 = F_zu / (pi * speicher_param["r_i_Gleitrohr"]**2)

    h_0 = h_WS - speicher_param["H_WS_OK_Dif"] - speicher_param["H_RS_Dif"]
    h_2 = h_0 - speicher_param["L_Fuehrung"]

    A_r = (pi * speicher_param["r_i_Fuehrung"]**2) - (pi * speicher_param["r_a_Doppelwand_Gleitrohr"]**2)
    A_b = (pi * speicher_param["r_Blende"]**2) - (pi * speicher_param["r_a_Doppelwand_Gleitrohr"]**2)
    A_1 = (pi * speicher_param["r_i_Gleitrohr"]**2)
    A_m = (pi * speicher_param["r_i_Fuehrung"]**2)
    A_1_stern = A_m - (pi * speicher_param["r_a_Einwand_Gleitrohr"]**2)
    A_0_stern = pi * 2 * speicher_param["R_Dif"] * speicher_param["H_RS_Dif"]

    i_0 = __Modell_find_index_h_pos(h_0, all_h_pos)
    i_2 = __Modell_find_index_h_pos(h_2, all_h_pos)
    i_1 = __Modell_find_index_h_pos(speicher_param["h_Rohrende"], all_h_pos)

    h_0_2 = all_h_pos[i_0] - all_h_pos[i_2]
    h_1_2 = all_h_pos[i_1] - all_h_pos[i_2]
    h_0_1 = all_h_pos[i_0] - all_h_pos[i_1]

    rho_0_2 = 0
    dh_0_2 = 0
    for i in range(i_2, i_0+1):
        hPos = all_h_pos[i]
        rho_0_2 += __Modell_Stoffwerte("rho",Speicherzustand[hPos][0]) * Speicherzustand[hPos][1]
        dh_0_2 += Speicherzustand[hPos][1]
    rho_0_2 /= dh_0_2

    theta_1_2 = 0
    dh_1_2 = 0
    for i in range(i_2, i_1+1):
        hPos = all_h_pos[i]
        theta_1_2 += Speicherzustand[hPos][0] * Speicherzustand[hPos][1]
        dh_1_2 += Speicherzustand[hPos][1]
    theta_1_2 /= dh_1_2

    theta_2 = Speicherzustand[all_h_pos[i_2]][0]
    f_WUE = 0.05
    theta_1_stern = (1 - f_WUE) * theta_2 + f_WUE * (Speicherzustand[all_h_pos[i_0]][0] + theta_1_2) / 2

    rho_2 = __Modell_Stoffwerte("rho", theta_2)
    rho_r = __Modell_Stoffwerte("rho", (theta_2 + theta_1_stern)/2)
    rho_1_stern = __Modell_Stoffwerte("rho", theta_1_stern)
    rho_w = rho_1

    zeta_B = 2.6
    zeta_M_0_stern = 1
    zeta_2_stern_1_stern = 1

    u_2_stern = u_1
    fak_neben_alt = 1
    wiederholen = True
    
    p_dyn = rho_1 * A_1 / A_m * u_1**2

    p_stat = 0
    alleTerme = 0
    
    while wiederholen:
        wiederholen = False
        p_stat = g * (rho_0_2 * h_0_2 - rho_r * h_1_2 - rho_w * h_0_1)
        u_1_stern = rho_2 * A_b * u_2_stern / (rho_1_stern * A_1_stern)
        u_M = (rho_1 * A_1 * u_1 + rho_1_stern * A_1_stern * u_1_stern) / (rho_w * A_m)
        u_0_stern = rho_w * A_m * u_M / (rho_w * A_0_stern)
        u_R = rho_2 * A_b * u_2_stern / (rho_r * A_r)

        alleTerme = p_stat + p_dyn - rho_w / 2 * u_M**2 + rho_1_stern * u_1_stern**2 * (A_1_stern/A_m - 0.5)\
                    - rho_w/2 * u_0_stern**2 * (1+zeta_M_0_stern) - zeta_2_stern_1_stern * rho_r/2 * u_R**2
        
        if alleTerme <= 0:
            fak_neben = 0
            break

        u_2_stern = (2/(rho_2*zeta_B)* alleTerme)**0.5
        fak_neben = u_2_stern / u_1 * A_b / A_1
        rho_w = (rho_1 + fak_neben *rho_2) / (fak_neben + 1)
        if abs(fak_neben_alt - fak_neben) > 1.0E-03:
            wiederholen = True
        fak_neben_alt = fak_neben
    m_neben_Plug_soll = dh_zu * fak_neben * rho_1
    m_neben_Plug = 0
    H_neben = 0
    counter_neben = 0
    while m_neben_Plug < m_neben_Plug_soll:
        rho_i_neben_Plug = __Modell_Stoffwerte("rho", Speicherzustand[all_h_pos[i_2+counter_neben]][0])
        m_i_neben_Plug = Speicherzustand[all_h_pos[i_2+counter_neben]][1] * rho_i_neben_Plug

        if m_neben_Plug_soll - m_neben_Plug <= m_i_neben_Plug:
            m_i_neben_Plug -= (m_neben_Plug_soll -m_neben_Plug)
            Speicherzustand[all_h_pos[i_2+counter_neben]][1] = m_i_neben_Plug / rho_i_neben_Plug
            H_neben += (m_neben_Plug_soll - m_neben_Plug) * __Modell_Stoffwerte("h", Speicherzustand[all_h_pos[i_2+counter_neben]][0])
            m_neben_Plug = m_neben_Plug_soll
        else:
            H_neben += m_i_neben_Plug * __Modell_Stoffwerte("h", Speicherzustand[all_h_pos[i_2+counter_neben]][0])
            m_neben_Plug += m_i_neben_Plug
            Speicherzustand[all_h_pos[i_2+counter_neben]][1] = 0
        counter_neben += 1
    H_neben += dh_zu * rho_1 * __Modell_Stoffwerte("h", theta_zu)
    m_zu_plus_neben = rho_1 * dh_zu * (1 + fak_neben)
    # NOTE: Hier kann eventuell etwas schlimmes passieren.
    theta_zu = __Modell_Stoffwerte("h_rev", H_neben/m_zu_plus_neben)
    rho_zu_plus_neben = __Modell_Stoffwerte("rho", theta_zu)
    dh_zu_plus_neben = m_zu_plus_neben / rho_zu_plus_neben
    dh_zu = dh_zu_plus_neben
    return dh_zu, theta_zu, fak_neben, Speicherzustand

# // Funktion: Nebenstrom ab
# \\ nur für Speicher Dessau notwendig
def __Modell_Nebenstrom_ab(F_ab, h_WS, dt, Speicherzustand):

    all_h_pos = sorted(list(Speicherzustand))
    F_neben = 0
    h_0 = h_WS - speicher_param["H_WS_OK_Dif"] - speicher_param["H_RS_Dif"]
    h_2 = h_0 - speicher_param["L_Fuehrung"]

    A_r = (pi * speicher_param["r_i_Fuehrung"]**2) - (pi * speicher_param["r_a_Doppelwand_Gleitrohr"]**2)
    A_b = (pi * speicher_param["r_Blende"]**2) - (pi * speicher_param["r_a_Doppelwand_Gleitrohr"]**2)
    A_1 = (pi * speicher_param["r_i_Gleitrohr"]**2)
    A_m = (pi * speicher_param["r_i_Fuehrung"]**2)
    A_1_stern = A_m - (pi * speicher_param["r_a_Einwand_Gleitrohr"]**2)
    A_0_stern = pi * 2 * speicher_param["R_Dif"] * speicher_param["H_RS_Dif"]

    i_0 = __Modell_find_index_h_pos(h_0, all_h_pos)
    i_2 = __Modell_find_index_h_pos(h_2, all_h_pos)
    i_1 = __Modell_find_index_h_pos(speicher_param["h_Rohrende"], all_h_pos)

    h_0_2 = all_h_pos[i_0] - all_h_pos[i_2]
    h_1_2 = all_h_pos[i_1] - all_h_pos[i_2]
    h_0_1 = all_h_pos[i_0] - all_h_pos[i_1]

    theta_2 = Speicherzustand[all_h_pos[i_2]][0]
    f_wue = 0.05 # fuer Dessau, 0.25 fuer Bautzen
    theta_1_stern = (1 - f_wue) * theta_2 + f_wue * Speicherzustand[all_h_pos[i_0]][0]

    rho_0_2 = 0
    dh_0_2 = 0
    for i in range(i_2, i_0+1):
        hPos = all_h_pos[i]
        rho_0_2 += __Modell_Stoffwerte("rho", Speicherzustand[hPos][0]) * Speicherzustand[hPos][1]
        dh_0_2 += Speicherzustand[hPos][1]
    rho_0_2 /= dh_0_2

    rho_2 = __Modell_Stoffwerte("rho", theta_2)
    rho_R = __Modell_Stoffwerte("rho", (theta_2 + theta_1_stern)/2)
    rho_1_stern = __Modell_Stoffwerte("rho", theta_1_stern)
    rho_w = __Modell_Stoffwerte("rho", Speicherzustand[all_h_pos[i_0]][0])

    zeta_B = 2.6 # fuer Dessau
    zeta_M_0_stern = 1
    zeta_2_stern_1_stern = 1

    u_1 = F_ab / A_1
    u_2_stern = u_1
    F_neben_alt = -1
    wiederholen = True

    p_dyn = rho_w * A_1 / A_m * u_1**2
    p_dyn += 0.75 # fuer Dessau

    counter = 0
    while wiederholen:
        counter += 1
        wiederholen = False
        p_stat = g * (rho_0_2 * h_0_2 - rho_R * h_1_2 - rho_w * h_0_1)
        u_1_stern = rho_2 * A_b * u_2_stern / (rho_1_stern * A_1_stern)
        u_R = rho_2 * A_b * u_2_stern / (rho_R * A_r)
        u_0_stern = (F_ab - F_neben) / A_0_stern
        u_M = (F_ab - F_neben) / A_m

        alleTerme = p_stat + p_dyn - rho_w/2 * u_M**2 - rho_1_stern * u_1_stern**2 *(A_1_stern/A_m+0.5)\
                    + rho_w/2 * u_0_stern**2 * zeta_M_0_stern - zeta_2_stern_1_stern * rho_R / 2 * u_R**2
        if alleTerme <= 0:
            F_neben = 0
            break
        u_2_stern = (2/(rho_2 * zeta_B)*alleTerme)**0.5
        F_neben = u_2_stern * A_b
        if abs(F_neben_alt - F_neben) > 1.0E-03:
            wiederholen = True
        F_neben_alt = F_neben
        
    F_neben = min(F_neben, 0.99 * F_ab)

    V_neben_soll = F_neben * dt
    V_neben_ist = 0
    m_neben = 0
    H_neben = 0
    counter_neben = 0
    while V_neben_ist < V_neben_soll:
        rho_i_Plug = __Modell_Stoffwerte("rho", Speicherzustand[all_h_pos[i_2+counter_neben]][0])
        h_i_Plug = __Modell_Stoffwerte("h", Speicherzustand[all_h_pos[i_2+counter_neben]][0])
        V_i_Plug = speicher_param["A_Quer"] * Speicherzustand[all_h_pos[i_2+counter_neben]][1]

        if V_neben_soll - V_neben_ist <= V_i_Plug:
            V_i_Plug -= (V_neben_soll - V_neben_ist)
            Speicherzustand[all_h_pos[i_2+counter_neben]][1] = V_i_Plug / speicher_param["A_Quer"]
            H_neben += (V_neben_soll - V_neben_ist) * rho_i_Plug * h_i_Plug
            m_neben += (V_neben_soll - V_neben_ist) * rho_i_Plug
            V_neben_ist += (V_neben_soll - V_neben_ist)
        else:
            H_neben += V_i_Plug * rho_i_Plug * h_i_Plug
            m_neben += V_i_Plug * rho_i_Plug
            Speicherzustand[all_h_pos[i_2+counter_neben]][1] = 0
            V_neben_ist += V_i_Plug

        counter_neben += 1
    
    theta_neben = 0
    if F_neben>0:
        theta_neben = __Modell_Stoffwerte("h_rev", H_neben/m_neben)
        rho_neben = __Modell_Stoffwerte("rho", theta_neben)
        F_neben = (m_neben/rho_neben) / dt #* 3600

    return theta_neben, F_neben, Speicherzustand

# // Funktion: find index of h_pos, muss nicht sein, aber ist halt so
def __Modell_find_index_h_pos(h_such, all_h_pos):
    ug = 0
    og = len(all_h_pos)
    while og - ug > 1:
        uneven = (ug + og) % 2
        i_mitte = (ug + og + uneven) / 2
        if h_such > all_h_pos[int(i_mitte)]:
            ug = int(i_mitte)
        else:
            og = int(i_mitte)
    if abs(h_such - all_h_pos[ug]) < abs(h_such - all_h_pos[og]):
        index = ug
    else:
        index = og
    return index

# // Funktion: Stoffwerte für bestimmte Temperatur ausgeben
def __Modell_Stoffwerte(groesse, theta=None):
    if groesse in ["rho","cp","lambda","TLF", "eta", "beta_rho", "h", "h_rev"]:
        Wert = __Temperatur_Abhaengige_Stoffwerte(groesse, theta)
    elif groesse in ["cp_Fundament", "rho_Fundament", "TLF_Fundament",
                     "lambda_Fundament", "rho_Mantel", "cp_Mantel"]:
        Wert = __Temperatur_Unabhaengige_Stoffwerte(groesse)
    else:
        raise ValueError("Die gesuchte Stoffgroesse '%s' ist nicht"%groesse
                         +"hinterlegt!\nAbbruch\n")
    return Wert

# // Funktion: Berechnung der Stoffwerte, die von der Temperatur abhängig sind
def __Temperatur_Abhaengige_Stoffwerte(groesse,theta):
    if groesse == "h_rev":
        h = theta / 1000
        Wert = -5.911685E-06 * h**2 + 2.420544E-01 * h - 4.700638E-01
        loop = 0
        while (loop < 2):
            h_theta_rev = __Modell_Stoffwerte("h", Wert)
            cp_theta_rev = __Modell_Stoffwerte("cp", Wert)
            Wert += (h * 1000 - h_theta_rev) / cp_theta_rev
            loop += 1
        return Wert

    # if (theta < 0) or (theta > 130):
    #     raise Warning("Temperaturwert von %s in __Modell_Stoffwerte() ausserhalb des zulaessigen Bereichs von 20 Grad C bis""130 Grad C. Gesuchte Groesse:" "'%s'.\n"%(theta,groesse))
    if groesse == "rho":
        return -2.525726E-03 * theta**2 - 2.123038E-01 * theta + 1.005011E+03
    elif groesse == "cp":
        return 9.776500E-03 * theta**2 - 7.677243E-01 * theta + 4.194836E+03
    elif groesse == "lambda":
        return (3.097195E-08 * theta**3 - 1.565775E-05 * theta**2 
                + 2.517120E-03 * theta + 5.531103E-01)
    elif groesse == "TLF":
        return -1.740136E-12 * theta**2 + 5.093712E-10 * theta+ 1.346697E-07
    elif groesse == "eta":
        return (-4.617641E-10 * theta**3 + 1.663679E-07 * theta**2 
                - 2.221812E-05 * theta + 1.301820E-03)
    elif groesse == "beta_rho":
        return 9.699776E-09 * theta**2 - 7.361887E-06 * theta - 1.135069E-04
    elif groesse == "h":
        return 4.394221E-01 * theta**2 + 4.129877E+03 * theta + 1.987100E+03

# // Funktion: Wert der Stoffwerte, die unabhängig von der Temperatur sind
def __Temperatur_Unabhaengige_Stoffwerte(groesse): 
    if groesse == "cp_Fundament":
        return 840
    elif groesse == "rho_Fundament":
        return 2400
    elif groesse == "lambda_Fundament":
        return 4
    elif groesse == "TLF_Fundament":
        return 2 / 840 / 2400
    elif groesse == "rho_Mantel":
        return 7800
    # J/kgK ?
    elif groesse == "cp_Mantel":
        return 490


def save_zustand(sz, name):
    os.makedirs("datei/sz/", exist_ok=True)
    with open("datei/sz/" + str(name)+".dat" , "w") as f:
        for key, v in sz.items():
            f.write(f"{float(key):.5f};{v[0]:.5f};{v[1]:.5f};{v[2]:.5f};{v[3]:.5f}\n")


