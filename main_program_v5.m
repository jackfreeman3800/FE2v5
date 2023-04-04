clear 
clc

rng default
load('LMM_model_calibrated');

students_decide_to_do_optional_additional_section = false;
num_simulations = 1;

for(k=1:num_simulations)


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Section:(Total 5 marks) You have to create a strategy to decide how
% much money is invested in the bonds (Assets) portfolio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%(For example -- LATAM)
t                               = 1;
AUM(t)                          = 4000000000; % AUM = Assets under management
fixed_amount_extract            = 0.01*AUM(t); % monthly fixed amount that investors extract every month
%%%========================= TO BE DONE! ======================
% Adjust the parameters using a strategy (5 Marks)
portion_to_maintain_in_cash     = 0.20;
%%%============================================================
actual_cash_amount              = AUM(t)*portion_to_maintain_in_cash;
actual_investment_amount        = AUM(t)*(1-portion_to_maintain_in_cash); % cash + bonds investment = AUM 
amount_to_subtract_interest_up  = 100000;
amount_to_add_interest_down     = 100000;



size_historical_term_structure = length(RateSpecG2); %From 2005 to 2013 


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Section:(Total 5 marks) You need to decide which periods are ran by 
% your simulation, and you have to "create" a strategy for
% including default risk (credit risk spread)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 3;      % num of years

%%%========================= TO BE DONE! ======================
% Random start selection -- You need to use a strategy for the selection of the
% first day between 10-May-2005 and 06-Jul-2011  (2 Marks)
current_day_ir_simulation = 1753; % 1 March 2010
%%%============================================================


list_days_simulation = busdays(RateSpecG2{1}.Settle+current_day_ir_simulation+1,RateSpecG2{1}.Settle+current_day_ir_simulation+365*N); % from next day to the next 3 years 


initial_day             = list_days_simulation(1);
current_month           = month(initial_day);
%%%========================= TO BE DONE! ======================
% The 10 year interest rate (r(0,10)) was selected for that random day (previous sub-section). Is that
% enough? HINT: What can be done to simulate Credit Spreads? (Remember the
% variable ``RateSpecG2'' contains the interest rate term structure of
% treasuries securities, that has different risk than the bonds that you
% might be using. (3 Marks)

R10Y_Note = getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + 10*365);


% Define a list of bond parameters
bond_prices = [131.28 133.11 83.43 92.09 101.34 113.53 105.85 109.17 123.44 118.93 96.02 104.97];
CouponRates = [0.085 0.0838 0.0775 0.06 0.0538 0.065 0.0445 0.053 0.08 0.0713 0.0555 0.065];
YearsToMaturity = [11 10 86 35 19 4 6 7 8 15 25 29];
fv = 1000; 

% Initialize empty lists to store bond yields and credit spreads
bond_yields = [];
credit_spreads = [];

% Iterate over the list of bond prices
for i = 1:length(bond_prices)
    % Get the bond data for the current iteration
    coupon_rate = CouponRates(i);
    years_to_maturity = YearsToMaturity(i);
    coupon_payment = fv * coupon_rate;
    total_payments = years_to_maturity * 2;
    price = bond_prices(i);
    
    % Compute the yield of the bond 
    bond_yield = (coupon_payment + (fv - price) / years_to_maturity) / ((fv + price) / 2);

    % Compute the credit spread as the difference between the bond yield and the treasury yield
    credit_spread = bond_yield - R10Y_Note;
    
    % Append the credit spread to its list
    credit_spreads = [credit_spreads credit_spread];
end

% Get an average of credit spreads to use for simulation
average_spread = mean(credit_spreads);

% Calculate an updated 'previous_rate_10Y_note'
previous_rate_10Y_note  = R10Y_Note + average_spread;

%%%============================================================

% (1/0.85)^(1/10)-1

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Section:(Total 10 marks) You have to adjust the suggested ``fictitious'' bonds with the ``real'' selected bonds, and 
% you have to decide the amount to be invested in every bond the first day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%========================= TO BE DONE! ======================
% Adjust the suggested ``fictitious'' bonds with the ``real'' selected bonds (TO BE DONE!) (10 Marks)
% For example, invest in two bonds with the following characteristics:
% Bond 1: KO 8 ½ 02/01/2022 Corp
% Bond 2: UPS 8 ⅜ 04/01/2020 Corp
% Bond 3: DOW 7 ¾ 10/01/2096 Corp
% Bond 4: GE 6 02/02/46 Corp
% Bond 5: DE 5 ⅜ 10/16/29 Corp
% Bond 6: CMCSA 6 ½ 01/15/15 Corp
% Bond 7: K 4.45 05/30/16 Corp
% Bond 8: MCD 5.3 03/15/17 Corp
% Bond 9: FDX 8 01/15/19 Corp
% Bond 10: LLY 7 ⅛ 06/01/25 Corp
% Bond 11: ALL 5.55 05/09/35 Corp
% Bond 12: MDLZ 6 ½ 02/09/40 Corp

P(1)          = 131.28; % HINT: The current prices of the bond reveal the credit spread with the treasuries rates. 
                    %       That difference can be used to simulate credit spreads along the simulation (The risk 
                    %       of the bonds that you use for the portfolio might be higher than treasuries.
P(2)          = 133.11;
P(3)          = 83.43;
P(4)          = 92.09;
P(5)          = 101.34;
P(6)          = 113.53;
P(7)          = 105.85;
P(8)          = 109.17;
P(9)          = 123.44;
P(10)          = 118.93;
P(11)          = 96.02;
P(12)          = 104.97;

coupon(1)     = 0.085;
coupon(2)     = 0.0838;
coupon(3)     = 0.0775;
coupon(4)     = 0.06;
coupon(5)     = 0.0538;
coupon(6)     = 0.065;
coupon(7)     = 0.0445;
coupon(8)     = 0.053;
coupon(9)     = 0.08;
coupon(10)     = 0.0713;
coupon(11)     = 0.0555;
coupon(12)     = 0.065;

face_value(1) = 1000;
face_value(2) = 1000;
face_value(3) = 1000;
face_value(4) = 1000;
face_value(5) = 1000;
face_value(6) = 1000;
face_value(7) = 1000;
face_value(8) = 1000;
face_value(9) = 1000;
face_value(10) = 1000;
face_value(11) = 1000;
face_value(12) = 1000;

bond_maturity(1) = 11;
bond_maturity(2) = 10;
bond_maturity(3) = 86;
bond_maturity(4) = 35;
bond_maturity(5) = 19;
bond_maturity(6) = 4;
bond_maturity(7) = 6;
bond_maturity(8) = 7;
bond_maturity(9) = 8;
bond_maturity(10) = 15;
bond_maturity(11) = 25;
bond_maturity(12) = 29;

bond_maturity_dates(1) = list_days_simulation(1)+bond_maturity(1)*365; % Make more accurate
bond_maturity_dates(2) = list_days_simulation(1)+bond_maturity(2)*365; 
bond_maturity_dates(3) = list_days_simulation(1)+bond_maturity(3)*365; 
bond_maturity_dates(4) = list_days_simulation(1)+bond_maturity(4)*365; % Make more accurate
bond_maturity_dates(5) = list_days_simulation(1)+bond_maturity(5)*365; 
bond_maturity_dates(6) = list_days_simulation(1)+bond_maturity(6)*365; 
bond_maturity_dates(7) = list_days_simulation(1)+bond_maturity(7)*365; % Make more accurate
bond_maturity_dates(8) = list_days_simulation(1)+bond_maturity(8)*365; 
bond_maturity_dates(9) = list_days_simulation(1)+bond_maturity(9)*365; 
bond_maturity_dates(10) = list_days_simulation(1)+bond_maturity(10)*365; % Make more accurate
bond_maturity_dates(11) = list_days_simulation(1)+bond_maturity(11)*365; 
bond_maturity_dates(12) = list_days_simulation(1)+bond_maturity(12)*365; 

coupons_per_year(1) = 2;
coupons_per_year(2) = 2;
coupons_per_year(3) = 2;
coupons_per_year(4) = 4;
coupons_per_year(5) = 2;
coupons_per_year(6) = 2;
coupons_per_year(7) = 2;
coupons_per_year(8) = 2;
coupons_per_year(9) = 2;
coupons_per_year(10) = 2;
coupons_per_year(11) = 2;
coupons_per_year(12) = 2;

bond_basis(1)       = 4; %Couple of these need to be changed
bond_basis(2)       = 4;
bond_basis(3)       = 4;
bond_basis(4)       = 4; 
bond_basis(5)       = 4;
bond_basis(6)       = 4;
bond_basis(7)       = 4; 
bond_basis(8)       = 4;
bond_basis(9)       = 4;
bond_basis(10)       = 4; 
bond_basis(11)       = 4;
bond_basis(12)       = 4;

w(1) = 1/12; % Portion invested in each bond
w(2) = 1/12; 
w(3) = 1/12; 
w(4) = 1/12; 
w(5) = 1/12; 
w(6) = 1/12; 
w(7) = 1/12; 
w(8) = 1/12; 
w(9) = 1/12; 
w(10) = 1/12; 
w(11) = 1/12; 
w(12) = 1/12; 

    
bond_units = [];
bond_units(t,:) = floor(w*actual_investment_amount./(P/100.*face_value)); % number of bonds to buy

amount_to_invest        = sum(bond_units(t,:).*P/100.*face_value);
actual_cash_amount      = actual_cash_amount + (actual_investment_amount-amount_to_invest); % We add the remaining cash from the rounding of the bond units to the cash account
cash(t)                 = actual_cash_amount;
actual_investment_amount= amount_to_invest; 
portfolio_bonds_current_value(t) = actual_investment_amount;
%%%============================================================

disp(['Day ' datestr(list_days_simulation(t)) ' -- AUM: $'  num2str(AUM(t)) ' -- Cash amount: ' num2str(cash(t)) ' -- Bonds Portfolio Value:' num2str(portfolio_bonds_current_value(t))]);

months_remaining_simulation  = 36; 

%%
t = 2; % We start the simulation from the second day


for(current_date = list_days_simulation(2:end)')
%for(current_date = list_days_simulation(1:30)')
    if(month(current_date)~= current_month)
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fourth Section: (Total 10 marks) Cash inflows/outflows (They happen at the
        % beginning of the month
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % We are starting a new month. We need to settle cash outflows (pension fund retirments) and
        % cash inflows (pension fund deposits), and update the value of the portfolio
        current_month = month(current_date);
        months_remaining_simulation = months_remaining_simulation-1; % we subtract one month from the remaining months of the simulation
        
        %%%========================= TO BE DONE! ======================
        % The 10 year interest rate (r(0,10)) was selected for that random day (previous sub-section). Is that
        % enough? (HINT: What can be done to simulate Credit Spreads? Remember the
        % variable ``RateSpecG2'' contains the interest rate term structure of
        % treasuries securities, that has different risk than the bonds that you
        % might be using. (4 Marks)

        R10Y_Note = getZeroRates(RateSpecG2{current_day_ir_simulation},RateSpecG2{current_day_ir_simulation}.Settle + 10*365);
        
        
        % Define a list of bond parameters
        bond_prices = [131.28 133.11 83.43 92.09 101.34 113.53 105.85 109.17 123.44 118.93 96.02 104.97];
        CouponRates = [0.085 0.0838 0.0775 0.06 0.0538 0.065 0.0445 0.053 0.08 0.0713 0.0555 0.065];
        YearsToMaturity = [11 10 86 35 19 4 6 7 8 15 25 29];
        fv = 1000; 
        
        % Initialize empty lists to store bond yields and credit spreads
        bond_yields = [];
        credit_spreads = [];
        
        % Iterate over the list of bond prices
        for i = 1:length(bond_prices)
            % Get the bond data for the current iteration
            coupon_rate = CouponRates(i);
            years_to_maturity = YearsToMaturity(i);
            coupon_payment = fv * coupon_rate;
            total_payments = years_to_maturity * 2;
            price = bond_prices(i);
            
            % Compute the yield of the bond 
            bond_yield = (coupon_payment + (fv - price) / years_to_maturity) / ((fv + price) / 2);
        
            % Compute the credit spread as the difference between the bond yield and the treasury yield
            credit_spread = bond_yield - R10Y_Note;
            
            % Append the credit spread to its list
            credit_spreads = [credit_spreads credit_spread];
        end
        
        % Get an average of credit spreads to use for simulation
        average_spread = mean(credit_spreads);
        
        % Calculate an updated 'previous_rate_10Y_note'
        current_rate_10Y_note  = R10Y_Note + average_spread;

        %%%============================================================
        
        if(current_rate_10Y_note-previous_rate_10Y_note>0.001) % Increase of of more than 10bp
            
            amount_to_subtract_interest_up  = ((current_rate_10Y_note-previous_rate_10Y_note)/0.001)*100000;
            
            if(actual_cash_amount-amount_to_subtract_interest_up>0) % We need to check we can extract such amount of cash
                
                
                
                actual_cash_amount = actual_cash_amount-(current_rate_10Y_note-previous_rate_10Y_note)*amount_to_subtract_interest_up;
            else
                %%%========================= TO BE DONE! ======================
                % You have to sell bonds to get cash for paying the
                % retirements of the fund (TO BE DONE!) (3 Marks)
                sales_proceeds = amount_to_subtract_interest_up - actual_cash_amount;    %NEW
                amount_to_invest = amount_to_invest - sales_proceeds;
                actual_cash_amount = actual_cash_amount + sales_proceeds;
                %%%============================================================
                
                % Then, you can subtract the cash
                actual_cash_amount = actual_cash_amount-amount_to_subtract_interest_up;
            end
                
        elseif(current_rate_10Y_note-previous_rate_10Y_note<-0.001) % Decrease of more than 10bp
            actual_cash_amount = actual_cash_amount+(previous_rate_10Y_note-current_rate_10Y_note)*amount_to_add_interest_down; % We don't need to check cash amount when it is added
             amount_to_add_interest_down     = ((previous_rate_10Y_note-current_rate_10Y_note)/0.001)*100000;

        else
            % Nothing changes in terms of the fixed cash inflows/outflows
            
        end
        if(actual_cash_amount-fixed_amount_extract>0)% We need to check we can extract such amount of cash
            actual_cash_amount = actual_cash_amount - fixed_amount_extract;
        else
            %%%========================= TO BE DONE! ======================
            % You have to sell bonds to get cash for paying the
            % retirements of the fund (TO BE DONE!) (3 Marks)
            sales_proceeds = fixed_amount_extract - actual_cash_amount;        %NEW
            amount_to_invest = amount_to_invest - sales_proceeds;
            actual_cash_amount = actual_cash_amount + sales_proceeds;
            %%%============================================================
            
            % Then, you can subtract the cash
            actual_cash_amount = actual_cash_amount-fixed_amount_extract; 
        end
        
        previous_rate_10Y_note  = current_rate_10Y_note;% update  current rate for the next step. Only changes once per month
       
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fifth Section (Total 60 marks): We change the composition of the portfolio: You
        % re-immunise the portfolio at the beginning of every month
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%========================= TO BE DONE! ======================        
        % Adjust the suggested ``fictitious'' bonds with the ``real'' selected bonds (TO BE DONE!) (5 Marks)
        Settle      = current_date;
        
        CouponRate  = coupon(1);        
        Maturity    = bond_maturity_dates(1);        
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(1);
        Period      = coupons_per_year(1);
        Basis       = bond_basis(1);
        [Price_t(1)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(1)        = mean(ZeroRates);
        
        CouponRate  = coupon(2);        
        Maturity    = bond_maturity_dates(2);        
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);        
        FaceValue   = face_value(2);
        Period      = coupons_per_year(2);
        Basis       = bond_basis(2);
        [Price_t(2)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(2)        = mean(ZeroRates);
        
        CouponRate  = coupon(3);
        Maturity    = bond_maturity_dates(3);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(3);
        Period      = coupons_per_year(3);
        Basis       = bond_basis(3);
        [Price_t(3)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(3)        = mean(ZeroRates);

        CouponRate  = coupon(4);
        Maturity    = bond_maturity_dates(4);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(4);
        Period      = coupons_per_year(4);
        Basis       = bond_basis(4);
        [Price_t(4)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(4)        = mean(ZeroRates);

        CouponRate  = coupon(5);
        Maturity    = bond_maturity_dates(5);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(5);
        Period      = coupons_per_year(5);
        Basis       = bond_basis(5);
        [Price_t(5)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(5)        = mean(ZeroRates);

        CouponRate  = coupon(6);
        Maturity    = bond_maturity_dates(6);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(6);
        Period      = coupons_per_year(6);
        Basis       = bond_basis(6);
        [Price_t(6)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(6)        = mean(ZeroRates);

        CouponRate  = coupon(7);
        Maturity    = bond_maturity_dates(7);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(7);
        Period      = coupons_per_year(7);
        Basis       = bond_basis(7);
        [Price_t(7)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(7)        = mean(ZeroRates);
        
        CouponRate  = coupon(8);
        Maturity    = bond_maturity_dates(8);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(8);
        Period      = coupons_per_year(8);
        Basis       = bond_basis(8);
        [Price_t(8)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(8)        = mean(ZeroRates);

        CouponRate  = coupon(9);
        Maturity    = bond_maturity_dates(9);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(9);
        Period      = coupons_per_year(9);
        Basis       = bond_basis(9);
        [Price_t(9)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(9)        = mean(ZeroRates);

        CouponRate  = coupon(10);
        Maturity    = bond_maturity_dates(10);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(10);
        Period      = coupons_per_year(10);
        Basis       = bond_basis(10);
        [Price_t(10)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(10)        = mean(ZeroRates);

        CouponRate  = coupon(11);
        Maturity    = bond_maturity_dates(11);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(11);
        Period      = coupons_per_year(11);
        Basis       = bond_basis(11);
        [Price_t(11)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(11)        = mean(ZeroRates);

        CouponRate  = coupon(12);
        Maturity    = bond_maturity_dates(12);
        ZeroDates   = [(current_date+1):365:Maturity];
        ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
        FaceValue   = face_value(12);
        Period      = coupons_per_year(12);
        Basis       = bond_basis(12);
        [Price_t(12)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
        r(12)        = mean(ZeroRates);

        %r_average = mean(r); %change to weighted average
        % Method 1 PV Liab
        % r_weighted_average = r(1)*w(1) + r(2)*w(2) + r(3)*w(3) + r(4)*w(4) + r(5)*w(5) + r(6)*w(6) + r(7)*w(7) + r(8)*w(8) + r(9)*w(9) + r(10)*w(10) +r(11)*w(11) +r(12)*w(12);
        % Method 2 - Iteration
        r_weighted_average =0;
        num_bonds = 12
        for i =1:num_bonds
            r_weighted_average = r_weighted_average + r(i)*w(i);
        end


        %%%============================================================
        
        %%%========================= TO BE DONE! ======================        
        % Improve the calculation of the duration (10 Marks)
        % coupon is list containing all bonds coupon rate,
        % Function used ([ModDuration,YearDuration,PerDuration] = bnddury(Yield,CouponRate,Settle,Maturity)
        % bnddury(Yield, CouponRate, Settle, Maturity, Period, Basis, EndMonthRule, 
        % IssueDate, FirstCouponDate, LastCouponDate, StartDate, Face)
        Settle = current_date;
        for(i=1:length(coupon))
            [ModDuration(i),Duration(i),PerDuration(i)] = bnddury(r(i),coupon(i),Settle,bond_maturity_dates(i),coupons_per_year(i),bond_basis(i),[],[],[],[],[],face_value(i));
        end
        %%%============================================================
        
        %%%========================= TO BE DONE! ======================        
        % Include convexity in the immunisation (25 Marks)
        %[YearConvexity,PerConvexity] = bndconvy(Yield,CouponRate,Settle,Maturity)
        for(i=1:length(coupon))
            [YearConvexity(i),PerConvexity(i)] = bndconvy(r(i),coupon(i),Settle,bond_maturity_dates(i),coupons_per_year(i),bond_basis(i),[],[],[],[],[],face_value(i));
        end
        
        num_total_bonds         = length(coupon);
        num_combination_bonds   = 6; % Here we select the number of bonds of our portfolio
        index_bond = nchoosek([1:num_total_bonds],num_combination_bonds); % This creates the list of ALL the possible combinations of bonds to immunise the portfolio
        
        %%%========================= TO BE DONE! ======================
        % Improve the estimation of the PV of Liability (8 Marks)        
        % months remaining
        % Calculate rolling weighted average of yields
        %Value_portfolio(j,:) = Price_t(index_bond(j,:)).*Units(j,:); % Might be in wrong place
        %weights = Value_portfolio(j,:) ./ sum(Value_portfolio(j,:)); % weights based on current portfolio value
        %weighted_r_average = sum(weights .* Yield(index_bond(j,:))); % weighted average of yields
        %Calculate PV of liabilty Using Weighted Average
        Liability           = months_remaining_simulation*fixed_amount_extract; % total fixed extracts per month remainings;
        Maturity_Portfolio  = months_remaining_simulation/12;
        m                   = 12; % compounds per year
        PV_Liability        =  Liability/(1+r_weighted_average/m).^(Maturity_Portfolio*m);
        
        
        for(j = 1:size(index_bond,1)) % WE test ALL the combinations. You have to select one from this for loop
            
            % Matrix for solutions
            
            A = [ones(1,num_combination_bonds); ...
                1/PV_Liability*[Duration(index_bond(j,:))], 
                1/PV_Liability*[YearConvexity(index_bond(j,:))]];
            
            b = [PV_Liability;...
                Maturity_Portfolio,
                Maturity_Portfolio^2];
            
            W(j,:) = A\b;
            
            Units(j,:) = round(W(j,:)./(Price_t(index_bond(j,:))));
            Value_portfolio(j,:) = Price_t(index_bond(j,:)).*Units(j,:);
            
            %Calculate Tracking Errror with convexity included
            %PV_Liability_Convexity= PV_Liability - 0.5 * sum (YearConvexity(index_bond(j,:)).* (Price_t(index_bond(j,:)) .* Units(j,:)).^2);
            TE(t,j) =  sum(Value_portfolio(j,:)) - PV_Liability;
        end
        %%%============================================================

        %%%========================= TO BE DONE! ======================        
        % In here you have to select the immunised portfolio.
        % We suggested the second, but this is not optimal.
        %%% EXTREMELY IMPORTANT: YOU HAVE TO EXPLAIN WHY YOU SELECTED SUCH
        %%% STRATEGY!!!!! (Include such explanation in the report)  (10 Marks)
        %%% ==========================================================
        immunisation_strategy_selected  = 3;
        if(~any(isnan(Units(immunisation_strategy_selected,:))))
            list_of_bonds_to_use            = index_bond(immunisation_strategy_selected,:);
            list_of_bonds_not_to_use        = setdiff(1:num_total_bonds,index_bond(immunisation_strategy_selected,:));
            bond_units(t,list_of_bonds_to_use) =  Units(immunisation_strategy_selected,:);
            bond_units(t,list_of_bonds_not_to_use) =  0;
            
            % Change in number of bonds
            cash_change = sum((bond_units(t,:)-bond_units(t-1,:)).*Price_t);
            actual_cash_amount = actual_cash_amount - cash_change;
        else
            % No change is made. It's not possible
            bond_units(t,:) = bond_units(t-1,:);
        end
        %%%============================================================

        if(students_decide_to_do_optional_additional_section)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Seventh Section: (Total 10 Marks): The current bond prices
            % provide a precise scenario of the interest rate term
            % structure. Nevertheless, the interest rate term structure is
            % dynamic. You might use the code below to improve the
            % immunisation code by including Monte Carlo simulations of the
            % possible future scenarios.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%========================= OPTIONAL! ======================
            
            nPeriods = 10;
            nTrials  = 10000;
            [LMMZeroRatesSimPaths, LMMForwardRatesSimPaths] = LMM{current_day_ir_simulation}.simTermStructs(nPeriods+1,'nTrials',nTrials,'antithetic',true);
            %%%============================================================
        end
    else
        %% ? (TO BE DONE?)
        
        
        bond_units(t,:) = bond_units(t-1,:);
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sixth Section: (Total 10 Marks): Updates in the value of the portfolio 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%========================= TO BE DONE! ======================
    % Adjust with the selected bonds (TO BE DONE!) (2 Marks)
    Settle      = current_date;
    
    CouponRate  = coupon(1);
    Maturity    = bond_maturity_dates(1);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);    
    FaceValue   = face_value(1);
    Period      = coupons_per_year(1);
    Basis       = bond_basis(1);
    [Price_t(1)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(1)        = mean(ZeroRates);

    CouponRate  = coupon(2);
    Maturity    = bond_maturity_dates(2);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(2);
    Period      = coupons_per_year(2);
    Basis       = bond_basis(2);
    [Price_t(2)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(2)        = mean(ZeroRates);
    
    CouponRate  = coupon(3);
    Maturity    = bond_maturity_dates(3);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);    
    FaceValue   = face_value(3);
    Period      = coupons_per_year(3);
    Basis       = bond_basis(3);
    [Price_t(3)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(3)        = mean(ZeroRates);

    CouponRate  = coupon(4);
    Maturity    = bond_maturity_dates(4);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(4);
    Period      = coupons_per_year(4);
    Basis       = bond_basis(4);
    [Price_t(4)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(4)        = mean(ZeroRates);

    CouponRate  = coupon(5);
    Maturity    = bond_maturity_dates(5);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(5);
    Period      = coupons_per_year(5);
    Basis       = bond_basis(5);
    [Price_t(5)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(5)        = mean(ZeroRates);

    CouponRate  = coupon(6);
    Maturity    = bond_maturity_dates(6);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(6);
    Period      = coupons_per_year(6);
    Basis       = bond_basis(6);
    [Price_t(6)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(6)        = mean(ZeroRates);

    CouponRate  = coupon(7);
    Maturity    = bond_maturity_dates(7);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(7);
    Period      = coupons_per_year(7);
    Basis       = bond_basis(7);
    [Price_t(7)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(7)        = mean(ZeroRates);
    
    CouponRate  = coupon(8);
    Maturity    = bond_maturity_dates(8);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(8);
    Period      = coupons_per_year(8);
    Basis       = bond_basis(8);
    [Price_t(8)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(8)        = mean(ZeroRates);

    CouponRate  = coupon(9);
    Maturity    = bond_maturity_dates(9);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(9);
    Period      = coupons_per_year(9);
    Basis       = bond_basis(9);
    [Price_t(9)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(9)        = mean(ZeroRates);

    CouponRate  = coupon(10);
    Maturity    = bond_maturity_dates(10);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(10);
    Period      = coupons_per_year(10);
    Basis       = bond_basis(10);
    [Price_t(10)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(10)        = mean(ZeroRates);

    CouponRate  = coupon(11);
    Maturity    = bond_maturity_dates(11);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(11);
    Period      = coupons_per_year(11);
    Basis       = bond_basis(11);
    [Price_t(11)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(11)        = mean(ZeroRates);

    CouponRate  = coupon(12);
    Maturity    = bond_maturity_dates(12);
    ZeroDates   = [(current_date+1):365:Maturity];
    ZeroRates   = RateSpecG2{current_day_ir_simulation}.getZeroRates(ZeroDates);
    FaceValue   = face_value(12);
    Period      = coupons_per_year(12);
    Basis       = bond_basis(12);
    [Price_t(12)] = prbyzero([Maturity CouponRate FaceValue Period Basis],Settle,ZeroRates,ZeroDates);
    r(12)        = mean(ZeroRates);
    
    %r_average = mean(r);
    % r_weighted_average = r(1)*w(1) + r(2)*w(2) + r(3)*w(3) + r(4)*w(4) + r(5)*w(5) + r(6)*w(6) + r(7)*w(7) + r(8)*w(8) + r(9)*w(9) + r(10)*w(10) +r(11)*w(11) +r(12)*w(12);
    % Method 2 - Iteration
    r_weighted_average =0;
        num_bonds = 12
        for i =1:num_bonds
            r_weighted_average = r_weighted_average + r(i)*w(i);
        end
    %%%============================================================
    
    % Current Valuation:
    portfolio_bonds_current_value(t)= sum(bond_units(t,:).*Price_t); % Where P_t is the price of bonds today (TO BE DONE!)
    cash(t)                         = actual_cash_amount;


    %%%========================= TO BE DONE! ======================
    % Improve the estimation of the PV of Liability (8 Marks)
    Liability                       = months_remaining_simulation*fixed_amount_extract; % total fixed extracts per month remainings;
    Maturity_Portfolio              = months_remaining_simulation/12;
    m                               = 12; % compounds per year        
    PV_Liability                    =  Liability/(1+r_weighted_average/m).^(Maturity_Portfolio*m);   
    %improve calculation by using the proportion invested in each bond
    % to calculate r_average rather than just getting the mean of the
    % interest rates
    %%%============================================================
    TE(t,1)                         =  portfolio_bonds_current_value(t) - PV_Liability;
    
    
    AUM(t) = portfolio_bonds_current_value(t) + cash(t);
    
    disp(['Day ' datestr(list_days_simulation(t)) ' -- AUM: $'  num2str(AUM(t)) ' -- Cash amount: ' num2str(cash(t)) ' -- Bonds Portfolio Value:' num2str(portfolio_bonds_current_value(t)) ' -- Tracking Error:' num2str(TE(t,1))]);

    % We increase the simulation day counters  by 1
    current_day_ir_simulation = current_day_ir_simulation+1;
    t = t+1;
end


TE_simulation{k} = TE;

end


