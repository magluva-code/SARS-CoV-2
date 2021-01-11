# Class Person

class Person(object):
    """
    An instance of Person is a part of the total population size N.
    Person p can have status as susceptible, Infected, Deceased, Recovered or Vaccinated.
    """
    def __init__(self, status="S"):
        # Needs custom exception.
        self.STATES = ["S", "I", "D", "R", "V"]
        if status not in self.STATES:
            raise Exception
        self.status = status
        self.symptomatic = False
        self.infectious = False
        self.quarantined = False
        self.history = {0 : self.status}

    def set_status(self, status, t):
        """
        @param: status (str), t (int)
        @return: void
        Change self.status and update self.history.
        """
        # Needs custom exception.
        if status not in self.STATES:
            raise Exception
        self.status = status
        self.update_history(t)

    def update_history(self, t):
        """
        @param: t (int)
        @return: void
        Update self.history.
        """
        self.history[t] = self.status

    def __call__(self):
        return self.status
