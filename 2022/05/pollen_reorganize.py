import pandas as pd


# Load in the Excel file as a DataFrame object.
path = "/Users/zamperini/Downloads/Flower pollen tracker data entry.xlsx"
df = pd.read_excel(path)

# Print out some information about the DataFrame.
print("Number of entries: {}".format(len(df)))
print("First 10 entries:\n{}".format(df.head(10)))

# This will be the new DataFrame where each entry in Part ID has been separated.
# Our new DataFrame will have the same columns as the old one.
new_df = pd.DataFrame(columns=df.columns)

def update_df(new_df, new_row):

    # Clean up some of the strings in each row.
    try:
        new_row["Event"] = new_row["Event"].lower()
    except:
        pass
    try:
        new_row["AM/PM"] = new_row["AM/PM"].lower()
    except:
        pass
    try:
        new_row["Crew "] = new_row["Crew "].upper()
    except:
        pass

    return pd.concat([new_df, new_row.to_frame().T])

# Calling iterrows() on a DataFrame in a loop like this mean each iteration of
# the loop is accessing the index (which is essentially the row number, starting
# at 0) and the row (which in pandas is called a "Series" object).
print("Creating new DataFrame...")
for index, row in df.iterrows():
    #print(index)

    # From inspection, the Excel file apparently has a bunch of blank rows
    # starting at index 168 and ending at index 195. Not a problem, just
    # skip those rows when we hit them with the "continue" command, which says
    # go to the next iteration of the loop.
    if index >= 168 and index <= 195:
        continue

    # When Event is "scape appears", there is no entry for Part ID. Since there
    # is no entry, using "split" below will cause an error. So handle those
    # cases up front by appending the row as-is to our DataFrame and then
    # moving onto the next row.
    if row["Event"] == "scape appears":
        new_df = update_df(new_df, row)
        continue

    # A single instance special case is row 584. It appears no data was taken.
    # We might as well just add it onto our new DataFrame anyways, but if you
    # don't want that then just comment out this operation (don't comment out
    # "continue" though).
    if index == 584:
        new_df = update_df(new_df, row)
        continue

    # Access the Part IDs from this Series, which is a string object.
    partids = row["Part ID"]

    # The split function will split a string into multiple strings, splitting
    # it at the specified character you tell it. It return a "list", which is
    # just a collection of each part ID each as an individual element.
    partids_split = partids.split(",")

    # We want to go through one partid at a time and create a new row (Series)
    # from it, keeping all else the same. We do this by copying the current
    # row we are on, and then overwriting Part ID with just the single ID.
    for partid in partids_split:
        new_row = row.copy()

        # A special case seems to be when mulitple consecutive Part IDs are
        # input, e.g., a-f. This is tricky, but still easily done if you know
        # what to Google. This is the trickiest part of the program, don't
        # stress over understanding this one.
        if len(partid.split("-")) > 1:
            start_letter, end_letter = partid.split("-")

            # Some special cases (y'all need to enter data more logically).
            # Ignore trying to understand this if you want.
            if index == 242:
                for chr_val in range(ord("a"), ord(end_letter)+1):
                    if chr_val == ord("a"):
                        new_row["Part ID"] = "aaa"
                        new_df = update_df(new_df, new_row)
                    else:
                        new_row["Part ID"] = chr(chr_val)
                        new_df = update_df(new_df, new_row)
            elif index == 505:
                for chr_val in range(ord("r"), ord("z")+1):
                    new_row["Part ID"] = chr(chr_val)
                    new_df = update_df(new_df, new_row)

                # The last two IDs for this entry.
                new_row["Part ID"] = "za"
                new_df = update_df(new_df, new_row)
                new_row["Part ID"] = "zb"
                new_df = update_df(new_df, new_row)

            # If entered with normal letters, e.g., a-g, this is how to handle it.
            else:
                for chr_val in range(ord(start_letter), ord(end_letter)+1):
                    new_row["Part ID"] = chr(chr_val)
                    new_df = update_df(new_df, new_row)

        # If not in these special cases, just overwrite with the single part ID.
        else:
            new_row["Part ID"] = partid
            new_df = update_df(new_df, new_row)

# The index (row numbers) are goofed up, so just reset them, then save to a new file.
new_df = new_df.reset_index()
new_df.to_excel("/Users/zamperini/Downloads/Flower pollen tracker data entry_edited.xlsx")
