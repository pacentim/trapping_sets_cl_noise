function subset = generate_error(n,wgt)
          % weight of error pattern
    last_generated_indicator=0;
    first_generated_indicator=0;
    if(last_generated_indicator==0)
        [subset,last_generated_indicator]=next_subset(n,wgt,first_generated_indicator);
        first_generated_indicator=last_generated_indicator;
    end
