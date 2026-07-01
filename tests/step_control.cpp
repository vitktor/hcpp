#include <gtest/gtest.h>
#include <hc/tracking/step_control.hpp>

using namespace hc;

TEST(AdaptiveStepController, GrowsAfterConsecutiveSuccesses) {
  // dt_init=0.1, growth kicks in after 3 consecutive successes, doubling dt
  AdaptiveStepController<double> controller(0.1, 0.001, 1.0);

  EXPECT_DOUBLE_EQ(controller.on_success(), 0.1); // streak 1, no growth yet
  EXPECT_DOUBLE_EQ(controller.on_success(), 0.1); // streak 2, no growth yet
  EXPECT_DOUBLE_EQ(controller.on_success(), 0.2); // streak 3, grows and resets
}

TEST(AdaptiveStepController, GrowthClampedAtMax) {
  AdaptiveStepController<double> controller(0.9, 0.001, 1.0);

  controller.on_success();
  controller.on_success();
  // raw growth would be 0.9 * 2 = 1.8, clamped to dt_max
  EXPECT_DOUBLE_EQ(controller.on_success(), 1.0);
}

TEST(AdaptiveStepController, ShrinksImmediatelyOnFailure) {
  // a single failure shrinks dt right away, no streak needed
  AdaptiveStepController<double> controller(0.1, 0.001, 1.0);

  EXPECT_DOUBLE_EQ(controller.on_failure(), 0.05);
}

TEST(AdaptiveStepController, ShrinkClampedAtMin) {
  AdaptiveStepController<double> controller(0.0015, 0.001, 1.0);

  // raw shrink would be 0.0015 * 0.5 = 0.00075, clamped to dt_min
  EXPECT_DOUBLE_EQ(controller.on_failure(), 0.001);
}

TEST(AdaptiveStepController, FailureResetsSuccessStreak) {
  AdaptiveStepController<double> controller(0.1, 0.0001, 10.0);

  controller.on_success(); // streak 1
  controller.on_success(); // streak 2
  EXPECT_DOUBLE_EQ(controller.on_failure(), 0.05); // streak reset to 0, dt halved

  // these two successes alone should NOT trigger growth, since the streak
  // needed to restart from zero after the failure
  controller.on_success(); // streak 1
  EXPECT_DOUBLE_EQ(controller.on_success(), 0.05); // streak 2, still no growth
}
