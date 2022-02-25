;; ####################
;; Pre-requisites
;; ####################


;;["reduce_transmissibility_to" [20 25 30]]
;;["agent_size" [0.02 0.03 0.04]]
;;["lockdown_from_day_n" 31]
;;["random_tests_percentage" 0 ]
;;["close_contact_tracing" true ]
;;["pcr_testing" true ]

extensions [nw csv] ;; csv extension used to read in in csv file

globals[
  exposedValues
  symptomaticValues
  severeValues
  criticalValues
  deathValues
  hours
  days
  SusceptibleAgents  ;; used for reporting purposes
  ExposedAgents  ;; used for reporting purposes
  InfectedAAgents  ;; used for reporting purposes
  InfectedMAgents  ;; used for reporting purposes
  InfectedSAgents  ;; used for reporting purposes
  InfectedCAgents  ;; used for reporting purposes
  RecoveredAgents  ;; used for reporting purposes
  DeadAgents  ;; used for reporting purposes
  isActiveTime?
  isRestTime?
  isWorkTime?
  isWeekend?
  isLevel5ActiveTime? ;; movement is only allowed during this time when in Level 5 lockdown
  PCR_sensitivity
  testing_capacity  ;; total number of tests per day (capacity)
  daily_CC_tests_remaining ;; remaining tests available today that can be used for close contact tracing (resets each day)
  daily_random_tests_percentage ;; percentage of testing budget assigned to random tests
  ireland_population
  ireland_testing_capacity
  in_lockdown?  ;; flag to indicate Level 5 lockdown is in place.  This will kick in from Day 30 onwards (i.e. December
  DaysList ;; list of days (day number)
  dailyConfirmedCasesList ;; List of confirmed cases per day
  dailyActualCasesList ;; List of actual new cases per day (number of people that became exposed per day)
  dailyActualCases ;; Counter for number of exposed people per day
  dailyConfirmedCases ;; Counter for number of confirmed cases per day
  dailyPCRTestsConducted ;; Counter for the number of PCR tests performed (excluding random tests)
  cumulativeConfirmedCases ;; Counter for cumulative number of confirmed cases
  cumulativeConfirmedRandomCases ;; counter of cumulative number of confirmed cases identified using random testing
  cumulativeNegativeRandomCases ;; counter of cumulative number of negative test results
  infectedPeople ;; agentset of all people infected (performance improvement)
  infectedPeopleNotIsolating ;; agentsert of all people infected but not isolating (performance improvement)
  SusceptiblePeople ;; agentset of all susceptible people (peformance improvement)
  casesValidation ;; used to store a list of real-world case data for model validation in behavioursearch tool
  model-adoption ;; used to store list of cases in model for model validation in behavioursearch tool
]

breed [households household] ;; household agents
breed [people person] ;; person agents

people-own[
  ;; each turtle has access to its own variables here
  epi_status ;; Susceptible Exposed Infected_A Infected_M Infected_S Infected_C Recovered Dead
  age
  gender
  isEssentialWorker?  ;; essential workers can continue to travel to work including during lockdown
  EconomicStatus  ;; employed, unemployed, student, retired
  person_id  ;; id of the person - can this be used as turtle id?
  family_id  ;; id of the famil (i.e. household agent) the person is a member of
  DaysExposed  ;; number of days agent is infected but not infectious
  DaysInfectiousToSymptomatic  ;; number of days agent is infectious
  DaysSevere  ;; number of days agent is severely ill in hospital
  DaysCritical ;;number of days agent is critically ill in ICU
  DaysCriticalToDeath  ;; number of days agent is critically ill before death
  DaysInfectiousToRecoveryAsymptomatic ;; number of days agent is infectious before recovery - asympotomatic cases
  DaysInfectiousToRecoveryMild  ;; number of days agent is infectious before recovery - mild cases
  DaysInfectiousToRecoverySevere  ;; number of days agent is infectious before recovery - severe cases
  DaysInfectiousToRecoveryCritical ;; number of days agent is infectious before recovery - critcal cases
  HomeLocation ;; coordinates of home location (agent)
  WorkLocation  ;; coordinates of work location (agent)
  days_since_infected ;;
  will_develop_symptoms?
  will_recover? ;; agent should only be checked once whether they will become severely / critically ill.  In the next iteration of the loop, they should not be checked again.
  will_goto_hospital?
  will_goto_icu?
  is_isolating? ;; agent will not infect anyone if they are isolating.  This flag is used to identify agents in quarantine
  close_contacts ;; list of contacts on a particular.  This list is reset each day
  tested_positive? ;; Has the person tested positive for COVID-19?
  is_close_contact? ;; Flags whether a person has been flagged as a close contact and will therefore be tested during the next testing iteration
  selected_for_random_test?
]

households-own[
  HouseholdType ;; one person, married or cohabiting couple / married or cohabiting couple with children / one parent family with children / two or more non-related persons
  Location ;; coordinates of household within patch
  isHospital?
  small_area_id
  small_area_name
  family_id
  isWorkplace?  ;; is a workplace where essential workers can go to work but is not a hospital
]

patches-own[
  patch_small_area ;; CSO small area assigned to this patch
]

;; ####################
;; Initialisation
;; ####################


;; This function takes mean and standard deviation parameters and returns a value from its log-normal distribution
;; mu = mean, sigma = standard deviation
to-report log-normal [mu sigma]
  let beta ln (1 + ((sigma ^ 2) / (mu ^ 2)))
  let x exp (random-normal (ln (mu) / (beta / 2)) sqrt beta)
  report mu
end

to setupGlobals
  ;; Disease dynamics age stratified probabilities
  set exposedValues [0.34 0.67 1 1 1 1 1 1.24 1.47 1.47]
  set symptomaticValues [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.9]
  set severeValues [0.0005 0.00165 0.00720 0.0208 0.0343 0.0765 0.1328 0.20655 0.2457 0.2457]
  set criticalValues [0.00003 0.00008 0.00036 0.00104 0.00216 0.00933 0.03639 0.08923 0.1742 0.1742]
  set deathValues [0.00002 0.00002 0.0001 0.00032 0.00098 0.00265 0.00766 0.02439 0.08292 0.16190]

  ;; initialise reporting variables
  set dailyConfirmedCasesList [0]
  set dailyActualCasesList [0]
  set dailyActualCases 0
  set dailyConfirmedCases 0
  set dailyPCRTestsConducted 0
  set daysList []
  set cumulativeConfirmedCases 0
  set cumulativeConfirmedRandomCases 0
  set cumulativeNegativeRandomCases 0
  set model-adoption []

  ;; initialise variables used in modeling
  set hours 0
  set days 0
  set PCR_sensitivity 0.87 ;; as per literature
  set ireland_population 4761865 ;; as per census 2016
  set ireland_testing_capacity 25000 ;; per day
  set testing_capacity ceiling ((ireland_testing_capacity / ireland_population) * count people) ;; set testing capacty to be % of the population (based on real world capacity of Ireland) - rounded up.
  set daily_random_tests_percentage random_tests_percentage
  set daily_CC_tests_remaining testing_capacity ;; reset daily testing capacity
  set in_Lockdown? False

  ;; Update some agentsets (performance improvement)
  set infectedPeople people with [ (epi_status = "Infected_A" or epi_status = "Infected_M" or epi_status = "Infected_S" or epi_status = "Infected_C") and is_isolating? = False]
  set infectedPeopleNotIsolating people with [epi_status = "Infected_A" or epi_status = "Infected_M"]
  set SusceptiblePeople people with [epi_status = "Susceptible"]

end


to import-smallareas
  let here_x min-pxcor
  let here_y max-pycor
  ;; read in file
  file-open "small_areas.txt"
  ;; loop through file and assign small area to patch
  while[not file-at-end?]
  [
    let items read-from-string (word "[" file-read-line "]")
    ask patch here_x here_y [ set patch_small_area item 0 items ]
    (ifelse
      here_x = max-pxcor and here_y = min-pycor [ stop ] ;; if we reach the last patch, quit

      here_x = max-pxcor [ ;; move onto next line
        set here_x min-pxcor
        set here_y here_y - 1
      ]
      ;; else
      [
        set here_x here_x + 1 ;; else, move one patch to the right
      ]
    )
  ]
 file-close
end

to import-households
  file-open "households.txt"

  while [not file-at-end?]
  [
    let items read-from-string (word "[" file-read-line "]")
    create-households 1 [
      ;; set variables from file
      set small_area_id item 1 items
      set small_area_name item 4 items
      set HouseholdType item 5 items
      set family_id item 6 items

      ;; set other variables
      set shape "house"
      set color pink
      set size Agent_size * 2
      set isHospital? false

      ;; move house to centre of correct patch, then set location within patch
      move-to one-of patches with [patch_small_area = [small_area_id] of myself]
      forward 0.1
      right random-normal 0 10
      left random-normal 1 10
      forward 0.1
    ]
  ]
  file-close

end

to import-people
  file-open "agents.txt"
  while [not file-at-end?]
  [
    let items read-from-string (word "[" file-read-line "]")

    create-people 1 [
      ;; set variables from file
      set gender item 2 items
      set age item 3 items
      set EconomicStatus item 4 items
      set isEssentialWorker? item 5 items
      set person_id item 6 items
      set family_id item 7 items

      ;; set other variables
      set shape "person"
      set color white
      set size Agent_size
      setxy random-xcor random-ycor
      set close_contacts []
      set tested_positive? False
      set is_close_contact? False
      set is_isolating? False
      ;; set all disease attributes for this agent
      InitialiseDiseaseAttributes

      ;; random-pxcor sets random coordinates to centre of a patch.  random-xcor sets a random set of coordinates
      setxy random-xcor random-ycor
    ]
  ]
  file-close
end


to InitialiseDiseaseAttributes
  set epi_status "Susceptible"
  set DaysExposed log-normal p_days_exposed 1.5 ;; default 4.5
  set DaysInfectiousToSymptomatic log-normal p_days_infectious_to_symptomatic 0.9 ;; default 1.1
  set DaysSevere log-normal p_days_severe 4.9 ;; default 6.6
  set DaysCritical log-normal p_days_critical 2  ;; default 1.5
  set DaysCriticalToDeath log-normal p_days_critical_to_death 4.8  ;; default 10.7
  set DaysInfectiousToRecoveryAsymptomatic log-normal p_days_infectious_to_recovery_asymptomatic 2  ;; default 8
  set DaysInfectiousToRecoveryMild log-normal p_days_infectious_to_recovery_mild 2  ;; default 8
  set DaysInfectiousToRecoverySevere log-normal p_days_infectious_to_recovery_severe 6.3 ;; default 18.1
  set DaysInfectiousToRecoveryCritical log-normal p_days_infectious_to_recovery_critical 6.3  ;;  default 18.1
  set days_since_infected 0

  ;; decide whether a person will develop symptoms or not if infected
  ifelse random-float 1 <= getProbability age symptomaticValues [
    ;; person will develop symptoms if infected
    set will_develop_symptoms? True
  ] ;; else
  [
    ;; else person will not
    set will_develop_symptoms? False
  ]

  ;; decide whether person will go to hospital if infected
  ifelse random-float 1 <= getProbability age severeValues [
    ;; person will become severely ill and go to hospital if infected
    set will_goto_hospital? True
  ] ;; else
  [
    ;; else person will not
    set will_goto_hospital? False
  ]

  ;; decide whether person will go to ICU if infected
  ifelse random-float 1 <= getProbability age criticalValues [
    ;; person will become critically ill and go to ICU if infected
    set will_goto_icu? True
  ] ;; else
  [
    ;; else person will not
    set will_goto_icu? False
  ]

  ;; decide whether person will recover from disease if infected
  ifelse random-float 1 <= getProbability age deathValues [
    ;; person will die if infected
    set will_recover? False
  ] ;; else
  [
    ;; else person will not
    set will_recover? True
  ]
end

to assign_people_to_households
  ask people[
    ;;set homeLocation one-of households with [isHospital? = false]
    set homeLocation one-of households with [family_id = [family_id] of myself]
  ]

  ;; Assign essential workers to a workplace
  ask people with [isEssentialWorker? = True] [
    set WorkLocation one-of households with [isWorkplace? = True]
  ]
end


to setup
  clear-all
  reset-ticks
  setupGlobals
  import-smallareas
  import-people
  import-households
  ;; Create hospital
  create-households 1[
    set color white
    set shape "hospital"
    set size Agent_size * 5
    set isHospital? true
    setxy random-xcor random-ycor
  ]
  ;; Create workplaces for essential workers
  create-households 100[
    set color blue
    set shape "factory"
    set size Agent_size * 5
    set isHospital? False
    set isWorkplace? True
    setxy random-xcor random-ycor
  ]
  assign_people_to_households
  ; infect a few random individuals
  ; n-of takes a sample of agents
  ask n-of 20 people[
    set epi_status "Exposed"
    set color grey
  ]
  setupGlobals
  read-in-data ;; read in sample data for validation only
end

;; ####################
;; Model Validation
;; ####################

;; used to load empirical data for validation
;; ref: https://www.youtube.com/watch?v=wMR-CEv1MVc
to read-in-data
  file-open "dailycases.csv"
  let header csv:from-row file-read-line
  set casesValidation []
  set model-adoption []
  while [ not file-at-end? ]
  [
    let data csv:from-row file-read-line
    let current first data
    if current > 0 [
      set casesValidation lput current casesValidation
    ]
  ]
  file-close
end


to-report compute-MSE [ series1 series2 ]
  ;; This function calculates the Mean Squared Error of modeled data vs real-world data
  ;; Only used by BehaviourSearch tool for model validation.  Not used in model itself
  ;; This is a measure of how different the time series are from each other.  We have the differences at each point, square them and then take the mean of these differences to get an average measure of the distance.
  ;; ref: https://www.youtube.com/watch?v=wMR-CEv1MVc
  let l1 length series1
  let l2 length series2
  let l min (list l1 l2)

  let minlist []
  let maxlist []

  if (l = l1) [ set minlist series1 set maxlist series2 ]
  if (l = l2) [ set minlist series2 set maxlist series1 ]

  set maxlist sublist maxlist 0 l

  let sqd-errors (map [ [?1 ?2] -> (abs ?1 - ?2) ^ 2 ] maxlist minlist)
  let mse (sum sqd-errors) / l

  report mse
end


;; ####################
;; Run the model
;; ####################

;; This reporter takes age and an epidemiology parameter list and returns the age
;; dependant probability for that parameter
to-report getProbability [#age #list]
  let x 0
  let itemNumber 0

  ;; convert age to a list index: Divide age by 10 and then convert to the largest integer less than that number
  ;; example 1: age 56 becomes 5.6 and then 5.  Index 5 in the list for symptomaticValues is 0.75
  ;; example 2: age 8 becomes 0.8 and then 0.  Index 0 in the list for symptomaticValues is 0.5
  set itemNumber #age / 10
  set itemNumber floor itemNumber
  set x item itemNumber #list

  report x
end

;; progress the disease depending on the time that has passed and the probability of progressing through each stage
to progress_epi_status
  ask people[
    (ifelse
      epi_status = "Susceptible" [
        ;; check age dependant probability
        set color white
      ]

      epi_status = "Exposed" [
        if days_since_infected >= DaysExposed[
          ;; check to see if agent becomes symptomatic
          ifelse will_develop_symptoms? = True [
            set epi_status "Infected_M"
            set color yellow
            set is_Isolating? True
          ][;; else agent is asymptomatic and will move into recovered compartment according to DaysInfectiousToRecoveryAsymptomatic
            set epi_status "Infected_A"
            set color brown
          ]
        ]
      ]
      ;; Progress from asymptomatic to recovered
      epi_status = "Infected_A"[
        if days_since_infected >= DaysExposed + DaysInfectiousToSymptomatic[
          if will_recover? = False [
            set epi_status "Infected_M"
            set color yellow
            set is_Isolating? True
          ]
          ;; Check if the disease has run its course
          if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic + DaysInfectiousToRecoveryAsymptomatic) [
            ;; agent recovers
            set epi_status "Recovered"
            set color green
            set close_contacts []
            set tested_positive? False
            set is_close_contact? False
            set is_isolating? False
          ]
        ]
      ]
      ;; Progress from Mild to Severe or recovered
      epi_status = "Infected_M"[
        if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic)[
          ;; check age dependant probability
          ifelse will_goto_hospital? = True or will_recover? = False [
            ;; agent gets severely ill and goes to hospital
            set epi_status "Infected_S"
            set color orange
            go_to_hospital self
          ][ ;; else agent will recover
            ;;set will_recover? true
            if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic + DaysInfectiousToRecoveryMild) [
              ;; agent recovers
              set epi_status "Recovered"
              set color green
              set close_contacts []
              set tested_positive? False
              set is_close_contact? False
              set is_isolating? False
            ]
          ]
        ]
      ]
      ;; Progress from Severe to Critical or recover
      epi_status = "Infected_S"[
        if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic + DaysSevere)[
          ;; check age dependant probability
          ifelse will_goto_icu? = True or will_recover? = False[
            ;; agent gets critically ill and goes to ICU
            set epi_status "Infected_C"
            set color red
          ] [;; else agent will recover
            ;; check this again - will a person die at this stage without becoming critical?
            if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic + DaysInfectiousToRecoverySevere)[
              ;; agent recovers
              set epi_status "Recovered"
              set color green
              set close_contacts []
              set tested_positive? False
              set is_close_contact? False
              set is_isolating? False
            ]
          ]
        ]
      ]

      epi_status = "Infected_C"[
        if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic + DaysSevere + DaysCritical)[
          ifelse will_recover? = false[
            ;; agent dies
            set epi_status "Dead"
            set color blue
          ][ ;; else agent recovers
            if days_since_infected >= (DaysExposed + DaysInfectiousToSymptomatic + DaysInfectiousToRecoveryCritical)[
              ;; agent recovers
              set epi_status "Recovered"
              set color green
              set close_contacts []
              set tested_positive? False
              set is_close_contact? False
              set is_isolating? False
            ]
          ]
        ]
      ]
    )
  ]
end

;; This procedure controls how people move around the environment depending on time of day and or what day of the week it is.
to MoveAgent
  ;; commute to/from work (for essential workers only)
  if isWeekend? = False and hours = 9 [go_to_work]
  if isWeekend? = False and hours = 17 [go_home_from_work]

  if isRestTime? = false[
    ;; 1. Behaviour of Retired or unemployed - they can roam freely during work hours
    ask people with [EconomicStatus = "Retired" or EconomicStatus = "Unemployed"] [
      if epi_status != "Dead" and epi_status != "Infected_S" and epi_status != "Infected_C" and is_isolating? = False[
        rt random-normal 0 10
        lt random-normal 0 10
        forward 0.4
      ]
    ]
    ;; 2. Behaviour of students
    ask people with [EconomicStatus = "Student"] [
      if epi_status != "Dead" and epi_status != "Infected_S" and epi_status != "Infected_C" and is_isolating? = False[
        rt random-normal 0 10
        lt random-normal 0 10
        forward 0.4
      ]
    ]
  ]

  ;; 3. Behaviour of working people - they don't move unless it's outside work hours and outside bed time
  if isWorkTime? = False and isRestTime? = False [
    ask people with [EconomicStatus = "Employed"] [
      if epi_status != "Dead" and epi_status != "Infected_S" and epi_status != "Infected_C" and is_isolating? = False[
        rt random-normal 0 10
        lt random-normal 0 10
        forward 0.4
      ]
    ]
  ]


  ;; Add schools / colleges and behaviour for students
  ;; Add behaviour for everyone during off time
end

;; This function controls agent behavious once in Level 5 lockdown
;; Including restriction of movement
to MoveAgentinLockdown
  ;; commute to/from work (for essential workers only)
  if isWeekend? = False and hours = 9 [go_to_work]
  if isWeekend? = False and hours = 17 [go_home]

  if isLevel5ActiveTime? = True [
    ask people with [epi_status != "Dead" and epi_status != "Infected_S" and epi_status != "Infected_C" and is_isolating? = False] [
      rt random-normal 0 10
      lt random-normal 0 10
      forward 0.4
    ]
  ]

end

;; This function performs a PCR test on a random selection of people regardless of whether or not they have been flagged as a close contact
;; People already symptomatic, dead or recovered are excluded
to runRandomTesting
  ;; calculate the number of random tests (rounded up) that can be performed
  let randomN ceiling (daily_random_tests_percentage * testing_capacity)
  ask n-of randomN people with [epi_status = "Susceptible" or epi_status = "Exposed" or epi_status = "Infected_A"] [
    ifelse random-float 1 <= PCR_sensitivity and epi_status != "Susceptible" [
      ;; test is positive
      set tested_positive? True
      set is_isolating? True
      set dailyConfirmedCases dailyConfirmedCases + 1
      set cumulativeConfirmedCases cumulativeConfirmedCases + 1
      set cumulativeConfirmedRandomCases cumulativeConfirmedRandomCases + 1
      ask other people with [member? who close_contacts = True][
        set is_close_contact? True
      ]
      set close_contacts []
      ] ;; else test is negative
      [
      set tested_positive? False
      set is_isolating? False
      set cumulativeNegativeRandomCases cumulativeNegativeRandomCases + 1
      ]
    set daily_CC_tests_remaining (daily_CC_tests_remaining - 1)
  ]

end

;; This function performs a PCR test on all close contacts, moves them to quarantine if required, and flags their close contacts for testing next time
;; Need to trigger value being set for is_close_contact or find some way to prompt a user
to runCCTesting
  ;; perform PCR test
  ;; If test is positive, do the following:
    ;; Send the person home to quarantine.  They cannot infect anyone while in quarantine.  Reset close contact list
    ;; Flag other people as close contact - they will quarantine and have PCR test which takes place when every 24 hours?
  ;; else
    ;;  send person out of quarantine
  ;; Perform the following for people either flagged as a close contact or if they have developed symptoms and are waiting for a test (i.e. are isolating already but have ot tested positive)
  ask people with [is_close_contact? = True or (epi_status = "Infected_M" and tested_positive? = False and is_isolating? = True)][
    ;; if testing capacity has not been reached, perform a PCR test
    if daily_CC_tests_remaining > 0 [ ;; perform the PCR test
      set dailyPCRTestsConducted dailyPCRTestsConducted + 1
      ifelse random-float 1 <= PCR_sensitivity and (epi_status = "Infected_A" or epi_status = "Infected_M" or epi_status = "Infected_S" or epi_status = "Infected_C") [
        ;; test is positive
        set tested_positive? True
        set is_isolating? True
        set dailyConfirmedCases dailyConfirmedCases + 1
        set cumulativeConfirmedCases cumulativeConfirmedCases + 1
        ask other people with [member? who close_contacts = True][
          set is_close_contact? True
        ]
        set close_contacts []
      ] ;; else test is negative
      [
        set tested_positive? False
        set is_isolating? False
      ]
      set daily_CC_tests_remaining (daily_CC_tests_remaining - 1)
    ]
  ]
end


;; Check when people get close to each other.  If they do, they may get infected
;; If people are already infected, they may have close contacts to report if they are confirmed positive
to check_close_contacts
  ; check for all infected if there are neighbours, and possibly infect them - or work the other way around - I will get infected from them
  ; count the number of infected turtles for each individual
  ; let creates a variable.  No = is used.  The variable name is the expression
  ; other turtles = all othere except for me
  ; with allows you to specify a condition
  ; can be restricted using in-radius (how close to me)
  ask people with [epi_status = "Susceptible"] [
    let infected-neighbours count other infectedPeople in-radius Agent_size
    if infected-neighbours > 0 [
      if random-float 1 <= getProbability age exposedValues and random 100 < reduce_transmissibility_to[
        set epi_status "Exposed"
        set color grey
        set days_since_infected 0
        set dailyActualCases dailyActualCases + 1
      ]
    ]
  ]

  ;; If Close contact tracing is enabled, track contacts
  if Close_contact_tracing [
    ask infectedPeopleNotIsolating [
      ;; if I am infected, record close contacts
      let susceptible-neighbours count other SusceptiblePeople in-radius Agent_size
      if susceptible-neighbours > 0 [
        set close_contacts lput [who] of other people in-radius Agent_size close_contacts
        set close_contacts reduce sentence close_contacts
        set close_contacts remove-duplicates close_contacts
      ]
  ]
  ]
end

to go
  if (all? people [epi_status != "Exposed" and epi_status != "Infected_A" and epi_status != "Infected_M" and epi_status != "Infected_S" and epi_status != "Infected_C"]) or days >= 96 [stop] ;; change days back to 96
  ifelse in_Lockdown? = True [
    MoveAgentinLockdown
  ]
  [ ;; else
    MoveAgent
  ]
  progress_epi_status

  check_close_contacts
  ;; Progress timing of disease for infected people
  ask people with [epi_status != "Susceptible" and epi_status != "Dead" and epi_status != "Recovered"] [
      ;;set days_since_infected days_since_infected + 0.083333 ;; 1 tick = 2 hours
      set days_since_infected days_since_infected + 0.041666 ;; 1 tick = 1 hour
    ]

  ;; Run random testing once per day.  The time is arbitary but happens before close contact testing
  if hours = 11 [runRandomTesting]

  ;; Run testing once per day to simulate a fixed amount of people that will get tested each day.  The time is arbitary
  if hours = 12 and PCR_Testing [runCCTesting]

  ;; Send all people that have tested positive or are isolating to go to quarantine (home)
  go_to_quarantine
  update_global_variables
end

;; This procedure sends the requested agent to a randomly selected hospital
;; The agent goes directly to the centre of that hospital agent.
to go_to_hospital [agentid]
  ask agentid [move-to one-of households with [isHospital? = true]]
end

;; This procedure sends the requested agent to its home
;; The agent goes directly to the centre of that building agent.
to go_home
  ask people with [epi_status != "Infected_S" and epi_status != "Infected_C"]  [move-to homeLocation]
end

;; This procedure sends the requested agent to its home
;; where they will isolate until further notice (14 days)
to go_to_quarantine
  ask people with [epi_status != "Infected_S" and epi_status != "Infected_C" and is_isolating? = True] [move-to homeLocation]
end

;; Procedure sends essential workers to their workplace
to go_to_work
  ask people with [isEssentialWorker? = True and is_Isolating? = False] [move-to WorkLocation]
end

to go_home_from_work
  ask people with [isEssentialWorker? = True and is_Isolating? = False] [move-to homeLocation]
end

to update_global_variables
  tick
  ;;set hours hours + 2
  set hours hours + 1
  ;;set days days + 0.083333
  set days days + 0.041666
  if hours >= 24 [
    set hours 0
    set dailyConfirmedCasesList lput (dailyConfirmedCases) dailyConfirmedCasesList ;; add number of confirmed cases to list.  There may be lots of other people infected but this is just confirmed cases
    set dailyActualCasesList lput dailyActualCases dailyActualCasesList
    set model-adoption lput dailyConfirmedCases model-adoption ;; used for model validation
    set DaysList lput round days DaysList
    set daily_CC_tests_remaining testing_capacity ;; reset daily testing capacity
    set dailyActualCases 0 ;; reset daily actual cases counter
    set dailyConfirmedCases 0 ;; resets daily conformed cases counter
    set dailyPCRTestsConducted 0 ;; resets daily PCR tests performed counter
  ]
  if days >= lockdown_from_day_n [
    set in_Lockdown? True
  ]



  ;; is today a weekend?
  ifelse (floor days mod 6 = 0) or (floor days mod 7 = 0) [
    set isWeekend? True
  ]
  [ ;; else
    set isWeekend? False
  ]

  ;; is this work time?
  ifelse hours >= 9 and hours <= 17 [
    set isWorkTime? True
  ]
  [ ;; else
    set isWorkTime? False
  ]

   ;; Is this active time?
  ifelse hours <= 7 or hours >= 21 [
    set isActiveTime? false
    set isRestTime? true
    if hours = 22 [go_home] ;; add if statement here to ensure this doesn't run every hour - resources!!
  ] [
    set isActiveTime? true
    set isRestTime? false
  ]

  ;; Is this active time in Level 5 lockdown?
  ifelse hours >= 18 and hours <= 21 [
    set isLevel5ActiveTime? True
  ]
  [
    set isLevel5ActiveTime? False
    if hours = 22 [go_home]
  ]

  ;; Update some agentsets (performance improvement
  set infectedPeople people with [ (epi_status = "Infected_A" or epi_status = "Infected_M" or epi_status = "Infected_S" or epi_status = "Infected_C") and is_isolating? = False]
  set infectedPeopleNotIsolating people with [epi_status = "Infected_A" or epi_status = "Infected_M"]
  set SusceptiblePeople people with [epi_status = "Susceptible"]


  ;; counter global variables used for results
  set SusceptibleAgents count people with [epi_status = "Susceptible"]
  set ExposedAgents count people with [epi_status = "Exposed"]
  set InfectedAAgents count people with [epi_status = "Infected_A"]
  set InfectedMAgents count people with [epi_status = "Infected_M"]
  set InfectedSAgents count people with [epi_status = "Infected_S"]
  set InfectedCAgents count people with [epi_status = "Infected_C"]
  set RecoveredAgents count people with [epi_status = "Recovered"]
  set DeadAgents count people with [epi_status = "Dead"]
end
@#$#@#$#@
GRAPHICS-WINDOW
9
10
667
869
-1
-1
50.0
1
10
1
1
1
0
1
1
1
-6
6
-8
8
1
1
1
ticks
30.0

BUTTON
667
10
756
78
Setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
756
10
848
79
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
665
80
1343
486
Full SEIR Model
Time (in hours)
Number of People
0.0
100.0
0.0
100.0
true
true
"" ""
PENS
"Susceptible" 1.0 0 -1513240 true "" "plot count people with [epi_status = \"Susceptible\"]"
"Exposed" 1.0 0 -7500403 true "" "plot count people with [epi_status = \"Exposed\"]"
"Infected Asymptomatic" 1.0 0 -6459832 true "" "plot count people with [epi_status = \"Infected_A\"]"
"Infected Mild" 1.0 0 -1184463 true "" "plot count people with [epi_status = \"Infected_M\"]"
"Infected Severe (Hospital)" 1.0 0 -955883 true "" "plot count people with [epi_status = \"Infected_S\"]"
"Infected Critical (ICU)" 1.0 0 -2674135 true "" "plot count people with [epi_status = \"Infected_C\"]"
"Dead" 1.0 0 -13345367 true "" "plot count people with [epi_status = \"Dead\"]"
"Recovered" 1.0 0 -13840069 true "" "plot count people with [epi_status = \"Recovered\"]"

MONITOR
1342
80
1399
125
Days
days
1
1
11

MONITOR
1401
80
1458
125
Hours
hours
1
1
11

SLIDER
847
10
1054
43
Reduce_transmissibility_to
Reduce_transmissibility_to
0
100
20.0
1
1
%
HORIZONTAL

SLIDER
847
43
1054
76
Agent_size
Agent_size
0
0.1
0.04
0.01
1
NIL
HORIZONTAL

MONITOR
1344
444
1468
489
Dead agents
DeadAgents
0
1
11

MONITOR
1342
127
1470
172
Susceptible agents
SusceptibleAgents
0
1
11

MONITOR
1342
173
1470
218
Exposed agents
ExposedAgents
0
1
11

MONITOR
1343
220
1470
265
Infected A agents
InfectedAAgents
0
1
11

MONITOR
1343
265
1469
310
Infected M agents
InfectedMAgents
0
1
11

MONITOR
1344
310
1469
355
Infected S agents
InfectedSAgents
0
1
11

MONITOR
1344
355
1469
400
Infected C agents
InfectedCAgents
0
1
11

MONITOR
1344
400
1469
445
Recovered agents
RecoveredAgents
0
1
11

PLOT
665
485
1346
861
Basic SEIR Model
Time (in hours)
Number of people
0.0
100.0
0.0
100.0
true
true
"" ""
PENS
"S - Susceptible" 1.0 0 -4539718 true "" "plot count people with [epi_status = \"Susceptible\"]"
"E - Exposed" 1.0 0 -11053225 true "" "plot count people with [epi_status = \"Exposed\"]"
"I - Infected" 1.0 0 -2674135 true "" "plot count people with [epi_status = \"Infected_M\"] + count people with [epi_status = \"Infected_A\"] + count people with [epi_status = \"Infected_S\"] + count people with [epi_status = \"Infected_C\"]"
"R - Removed" 1.0 0 -13840069 true "" "plot count people with [epi_status = \"Recovered\"] + count people with [epi_status = \"Dead\"]"

MONITOR
1748
487
1897
532
Daily Testing Capacity
testing_capacity
17
1
11

MONITOR
1751
579
1897
624
Tests remaining today
daily_CC_tests_remaining
0
1
11

MONITOR
1750
533
1909
578
People waiting for a test
count people with [is_close_contact? = True or (epi_status = \"Infected_M\" and tested_positive? = False and is_isolating? = True)]
0
1
11

MONITOR
1751
624
1894
669
Daily random tests
ceiling (daily_random_tests_percentage * testing_capacity)
17
1
11

MONITOR
1468
443
1655
488
Are we in Level 5 Lockdown?
in_Lockdown?
0
1
11

SLIDER
1054
10
1271
43
lockdown_from_day_n
lockdown_from_day_n
0
90
31.0
1
1
days
HORIZONTAL

PLOT
1346
668
1750
862
Daily Confirmed Cases
Time
Confirmed Cases
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot dailyConfirmedCases"

PLOT
1346
487
1750
670
Daily Actual Cases
Hours
Actual Cases
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot dailyActualCases"

SLIDER
1054
43
1270
76
random_tests_percentage
random_tests_percentage
0
1
0.5
0.1
1
NIL
HORIZONTAL

SWITCH
1271
10
1466
43
Close_contact_tracing
Close_contact_tracing
0
1
-1000

SWITCH
1272
44
1466
77
PCR_Testing
PCR_Testing
0
1
-1000

SLIDER
1940
51
2129
84
p_days_exposed
p_days_exposed
0
20
4.5
0.1
1
days
HORIZONTAL

SLIDER
1940
84
2243
117
p_days_infectious_to_symptomatic
p_days_infectious_to_symptomatic
0
20
1.1
0.1
1
days
HORIZONTAL

SLIDER
1940
118
2117
151
p_days_severe
p_days_severe
0
20
6.6
0.1
1
days
HORIZONTAL

SLIDER
1941
151
2119
184
p_days_critical
p_days_critical
0
20
1.5
0.1
1
days
HORIZONTAL

SLIDER
1942
184
2184
217
p_days_critical_to_death
p_days_critical_to_death
0
20
10.7
0.1
1
days
HORIZONTAL

SLIDER
1942
217
2313
250
p_days_infectious_to_recovery_asymptomatic
p_days_infectious_to_recovery_asymptomatic
0
20
8.0
0.1
1
days
HORIZONTAL

SLIDER
1942
249
2253
282
p_days_infectious_to_recovery_mild
p_days_infectious_to_recovery_mild
0
20
8.0
0.1
1
days
HORIZONTAL

SLIDER
1942
281
2267
314
p_days_infectious_to_recovery_severe
p_days_infectious_to_recovery_severe
0
20
18.1
0.1
1
days
HORIZONTAL

SLIDER
1942
314
2268
347
p_days_infectious_to_recovery_critical
p_days_infectious_to_recovery_critical
0
20
18.1
0.1
1
days
HORIZONTAL

TEXTBOX
1941
10
2233
60
Epidemiology Parameters
20
0.0
1

MONITOR
1751
669
1896
714
Positive random tests
cumulativeConfirmedRandomCases
0
1
11

MONITOR
1753
760
1902
805
Total confirmed cases
cumulativeConfirmedCases
0
1
11

MONITOR
1752
714
1903
759
Negative random tests
cumulativeNegativeRandomCases
0
1
11

MONITOR
1472
264
1629
309
Mild cases not isolating
;;count people with [epi_status = \"Infected_M\" and tested_positive? = False and is_isolating? = False]\ncount people with [epi_status = \"Infected_M\" and is_isolating? = False]
0
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

ambulance
false
0
Rectangle -7500403 true true 30 90 210 195
Polygon -7500403 true true 296 190 296 150 259 134 244 104 210 105 210 190
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Circle -16777216 true false 69 174 42
Rectangle -1 true false 288 158 297 173
Rectangle -1184463 true false 289 180 298 172
Rectangle -2674135 true false 29 151 298 158
Line -16777216 false 210 90 210 195
Rectangle -16777216 true false 83 116 128 133
Rectangle -16777216 true false 153 111 176 134
Line -7500403 true 165 105 165 135
Rectangle -7500403 true true 14 186 33 195
Line -13345367 false 45 135 75 120
Line -13345367 false 75 135 45 120
Line -13345367 false 60 112 60 142

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

building institution
false
0
Rectangle -7500403 true true 0 60 300 270
Rectangle -16777216 true false 130 196 168 256
Rectangle -16777216 false false 0 255 300 270
Polygon -7500403 true true 0 60 150 15 300 60
Polygon -16777216 false false 0 60 150 15 300 60
Circle -1 true false 135 26 30
Circle -16777216 false false 135 25 30
Rectangle -16777216 false false 0 60 300 75
Rectangle -16777216 false false 218 75 255 90
Rectangle -16777216 false false 218 240 255 255
Rectangle -16777216 false false 224 90 249 240
Rectangle -16777216 false false 45 75 82 90
Rectangle -16777216 false false 45 240 82 255
Rectangle -16777216 false false 51 90 76 240
Rectangle -16777216 false false 90 240 127 255
Rectangle -16777216 false false 90 75 127 90
Rectangle -16777216 false false 96 90 121 240
Rectangle -16777216 false false 179 90 204 240
Rectangle -16777216 false false 173 75 210 90
Rectangle -16777216 false false 173 240 210 255
Rectangle -16777216 false false 269 90 294 240
Rectangle -16777216 false false 263 75 300 90
Rectangle -16777216 false false 263 240 300 255
Rectangle -16777216 false false 0 240 37 255
Rectangle -16777216 false false 6 90 31 240
Rectangle -16777216 false false 0 75 37 90
Line -16777216 false 112 260 184 260
Line -16777216 false 105 265 196 265

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

factory
false
0
Rectangle -7500403 true true 76 194 285 270
Rectangle -7500403 true true 36 95 59 231
Rectangle -16777216 true false 90 210 270 240
Line -7500403 true 90 195 90 255
Line -7500403 true 120 195 120 255
Line -7500403 true 150 195 150 240
Line -7500403 true 180 195 180 255
Line -7500403 true 210 210 210 240
Line -7500403 true 240 210 240 240
Line -7500403 true 90 225 270 225
Circle -1 true false 37 73 32
Circle -1 true false 55 38 54
Circle -1 true false 96 21 42
Circle -1 true false 105 40 32
Circle -1 true false 129 19 42
Rectangle -7500403 true true 14 228 78 270

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

hospital
false
15
Rectangle -1 true true 75 45 105 255
Rectangle -1 true true 195 45 225 255
Rectangle -1 true true 105 135 195 165
Circle -1 false true 8 8 285

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="2. Cross validation experiment" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2500"/>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
  </experiment>
  <experiment name="0. test experiment" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count people with [epi_status = "Dead"]</metric>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="1. Iterations_experiment" repetitions="700" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
    <metric>cumulativeConfirmedCases</metric>
    <metric>cumulativeConfirmedRandomCases</metric>
    <metric>cumulativeNegativeRandomCases</metric>
    <enumeratedValueSet variable="Close_contact_tracing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reduce_transmissibility_to">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PCR_Testing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_exposed">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_from_day_n">
      <value value="31"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_symptomatic">
      <value value="1.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_severe">
      <value value="6.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_severe">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_asymptomatic">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Agent_size">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical_to_death">
      <value value="10.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_mild">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_critical">
      <value value="18.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="3. Real world data validation" repetitions="300" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
    <metric>cumulativeConfirmedCases</metric>
    <metric>cumulativeConfirmedRandomCases</metric>
    <metric>cumulativeNegativeRandomCases</metric>
    <enumeratedValueSet variable="Close_contact_tracing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reduce_transmissibility_to">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PCR_Testing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_exposed">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_from_day_n">
      <value value="31"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_symptomatic">
      <value value="1.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_severe">
      <value value="6.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_severe">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_asymptomatic">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Agent_size">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical_to_death">
      <value value="10.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_mild">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_critical">
      <value value="18.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="5. LIVE Experiment" repetitions="300" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
    <metric>cumulativeConfirmedCases</metric>
    <metric>cumulativeConfirmedRandomCases</metric>
    <metric>cumulativeNegativeRandomCases</metric>
    <metric>dailyPCRTestsConducted</metric>
    <enumeratedValueSet variable="Close_contact_tracing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reduce_transmissibility_to">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PCR_Testing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_exposed">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_from_day_n">
      <value value="31"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_symptomatic">
      <value value="1.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_severe">
      <value value="6.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_severe">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_asymptomatic">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Agent_size">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical_to_death">
      <value value="10.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_mild">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_critical">
      <value value="18.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="4.1 Sensitivity Analysis - lockdown" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
    <metric>cumulativeConfirmedCases</metric>
    <metric>cumulativeConfirmedRandomCases</metric>
    <metric>cumulativeNegativeRandomCases</metric>
    <metric>dailyPCRTestsConducted</metric>
    <enumeratedValueSet variable="Close_contact_tracing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reduce_transmissibility_to">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PCR_Testing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_exposed">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_from_day_n">
      <value value="10"/>
      <value value="30"/>
      <value value="50"/>
      <value value="70"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_symptomatic">
      <value value="1.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_severe">
      <value value="6.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_severe">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_asymptomatic">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Agent_size">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical_to_death">
      <value value="10.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_mild">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_critical">
      <value value="18.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="4.2 Sensitivity Analysis - testing and tracing" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
    <metric>cumulativeConfirmedCases</metric>
    <metric>cumulativeConfirmedRandomCases</metric>
    <metric>cumulativeNegativeRandomCases</metric>
    <metric>dailyPCRTestsConducted</metric>
    <enumeratedValueSet variable="Close_contact_tracing">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reduce_transmissibility_to">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PCR_Testing">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_exposed">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_from_day_n">
      <value value="31"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_symptomatic">
      <value value="1.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_severe">
      <value value="6.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_severe">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_asymptomatic">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Agent_size">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical_to_death">
      <value value="10.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_mild">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_critical">
      <value value="18.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="4.3 Sensitivity Analysis - reduce transmissibility" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>SusceptibleAgents</metric>
    <metric>ExposedAgents</metric>
    <metric>InfectedAAgents</metric>
    <metric>InfectedMAgents</metric>
    <metric>InfectedSAgents</metric>
    <metric>InfectedCAgents</metric>
    <metric>RecoveredAgents</metric>
    <metric>DeadAgents</metric>
    <metric>cumulativeConfirmedCases</metric>
    <metric>cumulativeConfirmedRandomCases</metric>
    <metric>cumulativeNegativeRandomCases</metric>
    <metric>dailyPCRTestsConducted</metric>
    <enumeratedValueSet variable="Close_contact_tracing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reduce_transmissibility_to">
      <value value="20"/>
      <value value="40"/>
      <value value="60"/>
      <value value="80"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PCR_Testing">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_exposed">
      <value value="4.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lockdown_from_day_n">
      <value value="31"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random_tests_percentage">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_symptomatic">
      <value value="1.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_severe">
      <value value="6.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_severe">
      <value value="18.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_asymptomatic">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Agent_size">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_critical_to_death">
      <value value="10.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_mild">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_days_infectious_to_recovery_critical">
      <value value="18.1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
