import hmmlearn
import hmmlearn.base



class phyloLL(hmmlearn.base._BaseHMM):
    def __init__(self):
        pass

    def _check(self):
        pass

    def _generate_sample_from_state(self, state, random_state=None):
        pass

    def _compute_log_likelihood(self, X):
        pass

    def _initialize_sufficient_statistics(self):
        pass

    def _accumulate_sufficient_statistics(self, stats, X, framelogprob,
                                          posteriors, fwdlattice, bwdlattice):
        pass

    def _do_mstep(self, stats):
        pass



# class Robot:
#     def __init__(self , name):
#         self.name = name
#     def say_hi(self):
#         print("Hi, I am " + self.name)
#     def say_bye(self):
#         print("see you!")
#
#
# class PhysicianRobot(Robot):
#     def say_hi(self):
#         print("Hello, this is " + self.name + ".  Nice to meet you!")
#     def say_bye(self):
#         print("it was good to see you, bye!")
#
#
# x = Robot("Mary")
# x.say_hi()
# x.say_bye()
#
# t = PhysicianRobot("Ali")
# t.say_hi()
# t.say_bye()